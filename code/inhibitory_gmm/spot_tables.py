"""Find a mouse's per-round mixed spot tables under /data and build the all-rounds cell x gene.

Three input layouts are accepted, tried in this order when ``source="auto"``:

``pairwise``      a pairwise-unmixing asset. Each ``<mouse>_R<N>/`` folder already holds the
                  mixed cell-by-gene tables that capsule built (``mixed_filtered_*`` and
                  ``mixed_all_spots_*``), so they are read as-is.
``spot_parquet``  an hcr-cache-spot-table asset: ``spots_R<N>.parquet`` + ``meta_R<N>.json``.
``processed``     one processed HCR asset per round, each with ``processing_manifest.json``
                  and ``image_spot_spectral_unmixing/mixed_spots_R<N>.pkl``.

Spot-level sources are turned into counts with the same rules as the pairwise capsule
(gene deletions, per-gene intensity floors, optional ``valid_spot`` gate), so every source
yields the same ``R<N>-<channel>-<gene>`` labelled pivot.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

SOURCES = ("auto", "pairwise", "spot_parquet", "processed")
SPOT_SUBSETS = ("filtered", "all_spots")

# Copied from hcr-pairwise-spot-unmixing so spot-level sources match its mixed tables.
SPOT_FILTER_GENES = {"Slc17a7": 200, "Slac17a7": 200}
GENE_DELETIONS = {"767022": {"R3": "Tac1"}, "755252": {"R3": "Tac1"}}
GENE_DICT_OVERRIDES = {
    "767022": {"R4": {"488": "Tac1", "514": "Crh", "561": "Calb1", "594": "Calb2", "638": "Npy"}},
}
MOUSE_ID_ALIASES = {"767018": ["767108"]}

_ROUND_DIR = re.compile(r"_R(\d+)$")
_LABEL = re.compile(r"^R(\d+)-(\w+?)-(.+)$")


@dataclass
class RoundInput:
    round: int
    kind: str  # "table" (pre-built long cell-by-gene) or "spots"
    path: Path
    source_asset: str
    genes: dict[str, str] = field(default_factory=dict)  # channel -> gene, spot sources only

    @property
    def key(self) -> str:
        return f"R{self.round}"


def _read_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def _mouse_ids(mouse_id: str) -> list[str]:
    return [mouse_id, *MOUSE_ID_ALIASES.get(mouse_id, [])]


def _gene_map(mouse_id: str, round_key: str, genes: dict[str, str]) -> dict[str, str]:
    genes = {str(c): g for c, g in genes.items()}
    genes.update(GENE_DICT_OVERRIDES.get(mouse_id, {}).get(round_key, {}))
    return genes


def find_pairwise_rounds(data_root: Path, mouse_id: str, spots: str) -> list[RoundInput]:
    assets = sorted(
        p for m in _mouse_ids(mouse_id) for p in data_root.glob(f"*{m}*_pairwise-unmixing*")
        if p.is_dir()
    )
    if not assets:
        return []
    asset = assets[-1]
    root = asset / "pairwise_unmixing" if (asset / "pairwise_unmixing").is_dir() else asset
    rounds = []
    for folder in sorted(p for p in root.iterdir() if p.is_dir()):
        match = _ROUND_DIR.search(folder.name)
        table = folder / f"mixed_{spots}_cell_by_gene.csv"
        if match and table.exists():
            rounds.append(RoundInput(int(match.group(1)), "table", table, asset.name))
    return sorted(rounds, key=lambda r: r.round)


def find_spot_parquet_rounds(data_root: Path, mouse_id: str) -> list[RoundInput]:
    rounds = []
    for asset in sorted(p for p in data_root.iterdir() if p.is_dir()):
        metas = sorted(asset.glob("meta_R*.json"))
        if not metas:
            continue
        # The asset name does not carry the mouse; subject.json (or data_description) does.
        record = next((p for p in (asset / "subject.json", asset / "data_description.json")
                       if p.exists()), None)
        if record and str(_read_json(record).get("subject_id")) not in _mouse_ids(mouse_id):
            continue
        for meta_path in metas:
            meta = _read_json(meta_path)
            parquet = asset / f"spots_{meta['key']}.parquet"
            if parquet.exists():
                rounds.append(RoundInput(int(meta["round"]), "spots", parquet, asset.name,
                                         _gene_map(mouse_id, meta["key"], meta.get("genes", {}))))
    return sorted(rounds, key=lambda r: r.round)


def find_processed_rounds(data_root: Path, mouse_id: str) -> list[RoundInput]:
    by_round: dict[int, RoundInput] = {}
    folders = sorted(
        p for m in _mouse_ids(mouse_id) for p in data_root.glob(f"HCR_{m}*_processed_*") if p.is_dir()
    )
    for folder in folders:
        manifest_path = next((p for p in (folder / "processing_manifest.json",
                                          folder / "derived" / "processing_manifest.json")
                              if p.exists()), None)
        if manifest_path is None:
            continue
        manifest = _read_json(manifest_path)
        if manifest.get("round") is None:
            continue
        n = int(manifest["round"])
        spots_dir = folder / "image_spot_spectral_unmixing"
        # Some datasets only have the stale R-1 name; the round always comes from the manifest.
        pkl = next((p for p in (spots_dir / f"mixed_spots_R{n}.pkl",
                                spots_dir / "mixed_spots_R-1.pkl") if p.exists()), None)
        if pkl is None:
            print(f"  {folder.name}: no mixed spot table, skipping")
            continue
        genes = {c: (v or {}).get("gene") for c, v in (manifest.get("gene_dict") or {}).items()}
        # Folders are name-sorted, so a later reprocess of the same round wins.
        by_round[n] = RoundInput(n, "spots", pkl, folder.name, _gene_map(mouse_id, f"R{n}", genes))
    return [by_round[n] for n in sorted(by_round)]


def discover_rounds(data_root: Path, mouse_id: str, spots: str = "all_spots",
                    source: str = "auto") -> tuple[str, list[RoundInput]]:
    """Return the source kind used and its rounds, trying each layout in order for ``auto``."""

    finders = {
        "pairwise": lambda: find_pairwise_rounds(data_root, mouse_id, spots),
        "spot_parquet": lambda: find_spot_parquet_rounds(data_root, mouse_id),
        "processed": lambda: find_processed_rounds(data_root, mouse_id),
    }
    for kind in (finders if source == "auto" else [source]):
        rounds = finders[kind]()
        if rounds:
            return kind, rounds
    raise FileNotFoundError(
        f"No mixed spot tables for mouse {mouse_id!r} (source={source!r}, spots={spots!r}) "
        f"under {data_root}"
    )


def _needed_columns(genes: dict[str, str]) -> list[str]:
    floors = [f"chan_{c}_intensity" for c, g in genes.items() if g in SPOT_FILTER_GENES]
    return ["chan", "cell_id", "valid_spot", *floors]


def _load_spots(path: Path, genes: dict[str, str]) -> pd.DataFrame:
    """Only the columns the counts need, with ``chan`` as a string categorical."""
    wanted = _needed_columns(genes)
    if path.suffix == ".parquet":
        import polars as pl

        frame = pl.read_parquet(path, columns=[c for c in wanted
                                               if c in pl.read_parquet_schema(path)])
        columns = {}
        for name in frame.columns:
            series = frame[name]
            if series.dtype == pl.Categorical:
                # Codes + categories avoid one Python string per spot.
                columns[name] = pd.Categorical.from_codes(
                    series.to_physical().to_numpy().astype("int32"),
                    categories=series.cat.get_categories().to_list(),
                )
            else:
                columns[name] = series.to_numpy()
        spots = pd.DataFrame(columns)
    else:
        spots = pd.read_pickle(path)
        spots = spots[[c for c in wanted if c in spots.columns]]
    spots["chan"] = spots["chan"].astype(str).astype("category")
    return spots


def spots_to_cell_by_gene(spots_df: pd.DataFrame, round_key: str, genes: dict[str, str],
                          mouse_id: str, spots: str) -> pd.DataFrame:
    """Long cell-by-gene (cell_id, gene, round_chan_gene, spot_count) from one round's spots."""

    if spots == "filtered":
        if "valid_spot" not in spots_df.columns:
            raise ValueError(
                f"{round_key}: spots='filtered' needs the valid_spot QC column, which only the "
                "pairwise-unmixing asset provides; use spots='all_spots' for this source."
            )
        spots_df = spots_df.loc[spots_df["valid_spot"].astype(bool)]
    chan = spots_df["chan"]
    if not isinstance(chan.dtype, pd.CategoricalDtype):
        chan = chan.astype(str).astype("category")
    chan = chan.cat.rename_categories(lambda c: str(c))
    spots_df = spots_df.assign(gene=chan.map(genes))
    spots_df = spots_df.dropna(subset=["gene"])

    deleted = GENE_DELETIONS.get(mouse_id, {}).get(round_key)
    if deleted:
        spots_df = spots_df.loc[spots_df["gene"] != deleted]

    channel_of = {g: c for c, g in genes.items() if g}
    for gene, floor in SPOT_FILTER_GENES.items():
        column = f"chan_{channel_of.get(gene)}_intensity"
        if gene in channel_of and column in spots_df.columns:
            spots_df = spots_df.loc[~((spots_df["gene"] == gene) & (spots_df[column] < floor))]

    counts = (
        spots_df.groupby(["cell_id", "gene"], observed=True)
        .size()
        .rename("spot_count")
        .reset_index()
    )
    counts["gene"] = counts["gene"].astype(str)
    counts.insert(2, "round_chan_gene",
                  round_key + "-" + counts["gene"].map(channel_of) + "-" + counts["gene"])
    return counts


def build_cell_by_gene(rounds: list[RoundInput], mouse_id: str, spots: str) -> pd.DataFrame:
    """All-rounds pivot: one row per cell, one ``R<N>-<channel>-<gene>`` column per gene."""

    pivots = []
    for r in rounds:
        print(f"  [{r.key}] {r.kind}: {r.path}")
        if r.kind == "table":
            long = pd.read_csv(r.path)
        else:
            long = spots_to_cell_by_gene(_load_spots(r.path, r.genes), r.key, r.genes, mouse_id, spots)
        pivots.append(long.pivot_table(index="cell_id", columns="round_chan_gene",
                                       values="spot_count", aggfunc="sum", fill_value=0))
    table = pd.concat(pivots, axis=1).fillna(0).astype(int)
    table.index.name = "cell_id"
    return table.rename(columns=lambda c: c.replace("Slac17a7", "Slc17a7"))


def sort_columns_by_round_channel(cxg: pd.DataFrame) -> pd.DataFrame:
    """Order ``R<N>-<channel>-<gene>`` columns by round, then channel; others go last."""

    def key(column: str):
        match = _LABEL.match(column)
        if not match:
            return (1, 0, 0, column)
        channel = match.group(2)
        return (0, int(match.group(1)), int(channel) if channel.isdigit() else 0, column)

    return cxg[sorted(cxg.columns, key=key)]


def plain_gene_table(cxg: pd.DataFrame) -> pd.DataFrame:
    """Collapse round-channel-gene columns to plain genes, summing genes imaged in >1 round."""

    return cxg.T.groupby(lambda c: c.rsplit("-", 1)[-1]).sum().T
