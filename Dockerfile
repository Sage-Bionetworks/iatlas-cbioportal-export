# uv + Python 3.10 preinstalled
FROM ghcr.io/astral-sh/uv:python3.10-bookworm
WORKDIR /root/cbioportal_export/

COPY . .

RUN uv sync --frozen --no-dev

# from here on, use uv run or call python in .venv explicitly
ENV PATH="/root/cbioportal_export/.venv/bin:$PATH"

WORKDIR /root/

# clone dep repos
RUN git clone https://github.com/rxu17/datahub-study-curation-tools.git -b upgrade-to-python3
RUN git clone https://github.com/cBioPortal/cbioportal.git -b v6.3.2

WORKDIR /root/cbioportal_export/