# hash:sha256:491b382563b5f4c88c54f3ff1a1a4658a3a2294e2ac49dd0b2a361f26027291a
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/code-server-python-extensions-pack:code-server4.108.2python3.12.4aind-metadata-mcp0.5.2withDOTFILES_REPOsupport

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-askpass /

RUN --mount=type=secret,id=git-access-token \
    GIT_ACCESS_TOKEN="$(cat /run/secrets/git-access-token)" \
    pip install -U --no-cache-dir \
    aind-data-access-api[docdb]==1.9.3 \
    codeocean==0.14.0 \
    plotly==6.6.0 \
    -e git+https://github.com/AllenInstitute/cell_type_mapper.git@79049fc679525ab79de2462ae23ed91ac2bfd554#egg=cell-type-mapper \
    -e git+https://github.com/AllenInstitute/xenium_analysis_tools@7844a9e5c3d2298316389e9392949db8fce0d34b#egg=xenium-analysis-tools \
    -e git+https://github.com/AllenNeuralDynamics/aind-hcr-data-loader.git@1cb418df118bf6833b552b4c05c3c809c375f2b7#egg=aind-hcr-data-loader \
    -e git+https://github.com/AllenNeuralDynamics/aind-hcr-qc.git@c7c3421b67fdfd555fdc7c3cddf852977b032407#egg=aind-hcr-qc
