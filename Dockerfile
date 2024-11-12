ARG MINIFORGE_VER="latest"

FROM condaforge/miniforge3:${MINIFORGE_VER}

WORKDIR /app

ADD . /app

RUN bash /app/shell/setup.sh

ENTRYPOINT ["/bin/bash", "/app/shell/AmpliPiper.sh"]
