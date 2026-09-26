"""Standard AIND metadata for the cell-typing derived asset (aind-data-schema 2.x).

Sources are the data assets the run actually read (the spot-parquet asset, a pairwise-unmixing
asset, or the processed rounds). The inheritance rules
(https://aind-data-schema.readthedocs.io/en/stable/inheritance.html) give two cases:

* one source whose name is ``<raw>_<process>_<ts>`` (e.g. a processed round): the asset is
  ``<raw>_cell-typing_<ts>`` and processing accumulates — ``Metadata.from_metadata``.
* several sources, or one source that already merges several acquisitions (the spot-parquet
  and pairwise-unmixing assets have no single raw input): the ANALYZED name
  ``cell-types-and-learning_cell-typing_<ts>``, subject kept (one mouse), instrument and
  acquisition dropped, processing started fresh.

Legacy v1 records are upgraded with aind-metadata-upgrader first; a round whose subject.json
cannot be upgraded (a placeholder with no species) contributes its data_description only.
"""
from __future__ import annotations

import json
import re
import warnings
from datetime import datetime, timezone
from pathlib import Path

from aind_data_schema.components.identifiers import Code
from aind_data_schema.core.data_description import DataDescription
from aind_data_schema.core.metadata import Metadata
from aind_data_schema.core.processing import DataProcess, Processing, ProcessStage
from aind_data_schema.core.subject import Subject
from aind_data_schema.utils.inheritance import derive_data_description_analyzed
from aind_data_schema_models.data_name_patterns import DataRegex, datetime_to_name_string
from aind_data_schema_models.process_names import ProcessName
from aind_metadata_upgrader.data_description.v1v2 import DataDescriptionV1V2
from aind_metadata_upgrader.upgrade import Upgrade

PROCESS_NAME = "cell-typing"
# Prefix for multi-source (ANALYZED) names; matches hcr-cache-spot-table's.
ASSET_NAME_PREFIX = "cell-types-and-learning"
CODE_URL = "https://github.com/AllenNeuralDynamics/hcr-unmix-to-cells"
CAPSULE_ID = "ae35f2e0-5603-4404-8ed4-a22c5a6fd421"
EXPERIMENTERS = ["Matt Davis"]
SOFTWARE_VERSION = "0.5.0"


def _raw_input(name: str) -> str | None:
    """The raw asset a ``<raw>_<process>_<ts>`` name derives from, or None if it has none."""
    match = re.match(DataRegex.DERIVED.value, name)
    return match.group("input") if match else None


def asset_name(source_names: list[str], creation_time: datetime) -> str:
    """The derived asset name the rules give for these sources (the launcher uses it too)."""
    stamp = datetime_to_name_string(creation_time)
    if len(source_names) == 1 and (raw := _raw_input(source_names[0])):
        return f"{raw}_{PROCESS_NAME}_{stamp}"
    return f"{ASSET_NAME_PREFIX}_{PROCESS_NAME}_{stamp}"


def _read(path: Path) -> dict | None:
    if not path.exists():
        return None
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _upgrade_v1(folder: Path, record: dict, log=print) -> Metadata:
    record["data_description"].pop("input_data_name", None)  # else resolved through DocDB
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return Upgrade(record).metadata
        except Exception as exc:  # e.g. a placeholder subject.json with no species name
            if "subject" not in record:
                raise
            log(f"  WARNING: {folder.name}/subject.json cannot be upgraded ({exc!r}); "
                "subject omitted, subject_id stays in data_description")
            target = DataDescription.model_fields["schema_version"].default
            dd = DataDescriptionV1V2().upgrade(record["data_description"], target,
                                               metadata={"name": folder.name})
            dd = DataDescription.model_validate(dd)
            return Metadata.model_construct(name=dd.name, location=str(folder),
                                            data_description=dd)


def load_source(folder: Path, log=print) -> Metadata:
    """One mounted source asset's metadata: data_description + subject (+ processing if v2)."""
    folder = Path(folder)
    dd = _read(folder / "data_description.json")
    if dd is None:
        raise SystemExit(f"ERROR: {folder.name} has no data_description.json to inherit from.")
    if str(dd.get("schema_version", "")).split(".")[0] in ("0", "1"):
        log(f"  {folder.name}: upgrading v{dd.get('schema_version')} metadata")
        record = {"name": folder.name, "location": str(folder), "data_description": dd}
        subject = _read(folder / "subject.json")
        if subject is not None:
            record["subject"] = subject
        return _upgrade_v1(folder, record, log)
    fields = {"data_description": DataDescription.model_validate(dd)}
    if (subject := _read(folder / "subject.json")) is not None:
        fields["subject"] = Subject.model_validate(subject)
    if (processing := _read(folder / "processing.json")) is not None:
        fields["processing"] = Processing.model_validate(processing)
    name = fields["data_description"].name
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return Metadata(name=name, location=str(folder), **fields)
        except Exception:  # no subject/processing/model: fine as a source, not as an asset
            return Metadata.model_construct(name=name, location=str(folder), **fields)


def build(source_folders, parameters: dict, start: datetime, outputs: dict, summary: str,
          subject_id: str, creation_time: datetime | None = None, log=print) -> Metadata:
    end = datetime.now(timezone.utc)
    creation_time = creation_time or end
    sources = [load_source(f, log) for f in source_folders]
    names = [m.data_description.name for m in sources]
    subject_ids = sorted({str(m.data_description.subject_id) for m in sources})
    if subject_ids != [str(subject_id)]:
        raise SystemExit(f"ERROR: sources are for subjects {subject_ids}, not {subject_id}.")
    name = asset_name(names, creation_time)
    single_raw = len(sources) == 1 and _raw_input(names[0]) is not None
    pattern = DataRegex.DERIVED if single_raw else DataRegex.ANALYZED
    if not re.match(pattern.value, name):  # pragma: no cover
        raise SystemExit(f"ERROR: asset name {name!r} does not match {pattern.name}.")

    processing = Processing(data_processes=[DataProcess(
        process_type=ProcessName.ANALYSIS,
        name="Cell typing (inhibitory GMM classes, subclasses, subtypes)",
        stage=ProcessStage.ANALYSIS,
        code=Code(url=CODE_URL, version=SOFTWARE_VERSION,
                  parameters={**parameters, "capsule_id": CAPSULE_ID}),
        experimenters=list(EXPERIMENTERS),
        start_date_time=start,
        end_date_time=end,
        output_path=".",
        output_parameters=outputs,
        notes="Procedure and column definitions: cell_typing_table.md.",
    )])
    dd_kwargs = dict(creation_time=creation_time, tags=["HCR", PROCESS_NAME, str(subject_id)],
                     data_summary=summary)
    location = f"s3://aind-open-data/{name}"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if single_raw:
            derived = Metadata.from_metadata(sources[0], process_name=PROCESS_NAME,
                                             location=location, new_processing=processing,
                                             **dd_kwargs)
            record = derived.model_dump(mode="json")
        else:
            # Multi-acquisition: one subject, so subject is inherited (from the first source
            # that has one); instrument / acquisition / procedures-per-round are dropped.
            dd = derive_data_description_analyzed(sources[0].data_description, PROCESS_NAME,
                                                  source_data=names, **dd_kwargs)
            subject = next((m.subject for m in sources if m.subject is not None), None)
            record = {"name": name, "location": location,
                      "data_description": dd.model_dump(mode="json"),
                      "subject": subject.model_dump(mode="json") if subject else None,
                      "processing": processing.model_dump(mode="json")}
        # Pin the name (our prefix, not the inherited project_name) and re-validate the record.
        record["name"] = record["data_description"]["name"] = name
        record["location"] = location
        derived = Metadata.model_validate(record)
    return derived


def write(results_dir, source_folders, parameters: dict, start: datetime, outputs: dict,
          summary: str, subject_id: str, creation_time: datetime | None = None,
          log=print) -> str:
    """Write the derived asset's core metadata files into *results_dir*; return its name."""
    derived = build(source_folders, parameters, start, outputs, summary, subject_id,
                    creation_time=creation_time, log=log)
    res = Path(results_dir)
    derived.write_standard_files(output_directory=res)
    derived.write_standard_file(output_directory=res)  # metadata.nd.json
    dd = derived.data_description
    log(f"[aind_metadata] {dd.name}  subject_id={dd.subject_id}  project={dd.project_name}")
    log(f"[aind_metadata] source_data={dd.source_data}  "
        f"processes={[p.name for p in derived.processing.data_processes]}")
    return dd.name
