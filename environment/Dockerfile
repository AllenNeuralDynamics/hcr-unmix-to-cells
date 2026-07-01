# hash:sha256:06a488d4776f77bb0f7d4b21c9126f254d2b8af0b67cffc97cb3cfdbd0c7644b
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/code-server-python-extensions-pack:code-server4.108.2python3.12.4aind-metadata-mcp0.5.2withDOTFILES_REPOsupport

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

RUN --mount=type=secret,id=git-access-token \
    GIT_ACCESS_TOKEN="$(cat /run/secrets/git-access-token)" \
    pip install -U --no-cache-dir \
    aind-data-access-api[docdb]==1.9.3 \
    anndata==0.12.18 \
    codeocean==0.14.0 \
    igraph==1.0.0 \
    leidenalg==0.12.0 \
    plotly==6.6.0 \
    polars==1.42.0 \
    scanpy==1.12.2 \
    seaborn==0.13.2 \
    -e git+https://github.com/AllenInstitute/cell_type_mapper.git@79049fc679525ab79de2462ae23ed91ac2bfd554#egg=cell-type-mapper \
    -e git+https://github.com/AllenInstitute/xenium_analysis_tools@7844a9e5c3d2298316389e9392949db8fce0d34b#egg=xenium-analysis-tools \
    -e git+https://github.com/AllenNeuralDynamics/aind-hcr-data-loader.git@0f93c97bba31c5d6ffdacd656f3cc159a9ceae9f#egg=aind-hcr-data-loader \
    -e git+https://github.com/AllenNeuralDynamics/aind-hcr-qc.git@e03778d64390b3e70c7cea81ac923944520a55eb#egg=aind-hcr-qc \
    -e git+https://github.com/AllenNeuralDynamics/ophys-mfish-dataset-catalog.git@d54316aa156fa8302e87615573a9d5cc89e58fa1#egg=ophys-mfish-dataset-catalog

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
