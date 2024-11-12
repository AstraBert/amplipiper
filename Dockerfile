ARG CONDA_VER=latest
ARG MAMBA_VER=latest


FROM mambaorg/micromamba:${MAMBA_VER} as mambabase
FROM conda/miniconda3:${CONDA_VER}

COPY --from=mambabase /usr/bin/micromamba /usr/bin/

RUN echo "deb http://archive.debian.org/debian stretch main contrib non-free" > /etc/apt/sources.list && apt update && apt install -y git

WORKDIR /app

ADD . /app

RUN bash /app/shell/setup.sh

ENTRYPOINT ["/bin/bash", "/app/shell/AmpliPiper.sh"]
