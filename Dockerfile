# uv + Python 3.10 preinstalled
FROM ghcr.io/astral-sh/uv:python3.10-bookworm
WORKDIR /workspace/iatlas_cbioportal_export

COPY . .

RUN uv sync --frozen --no-dev

# make sure any UID can read/traverse
RUN chmod -R a+rX /workspace/iatlas_cbioportal_export

# from here on, use uv run or call python in .venv explicitly
ENV PATH="/workspace/iatlas_cbioportal_export/.venv/bin:$PATH"

WORKDIR /workspace/

# clone dep repos
RUN git clone https://github.com/rxu17/datahub-study-curation-tools.git -b upgrade-to-python3
RUN git clone https://github.com/cBioPortal/cbioportal.git -b v6.3.2

RUN chmod -R a+rX /workspace/datahub-study-curation-tools /workspace/cbioportal


WORKDIR /workspace/iatlas_cbioportal_export