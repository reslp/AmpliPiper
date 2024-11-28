ARG MINIFORGE_VER="latest"

FROM condaforge/miniforge3:${MINIFORGE_VER}

WORKDIR /app

RUN git clone https://github.com/nhmvienna/AmpliPiper.git

RUN bash /app/AmpliPiper/shell/setup.sh

RUN chmod +x /app/AmpliPiper/shell/AmpliPiper.sh

ENV PATH=$PATH:/app/AmpliPiper/shell
RUN sed -i 's|tmp=$(dirname $0)|tmp=/app/AmpliPiper/shell|' /app/AmpliPiper/shell/AmpliPiper.sh

CMD ["AmpliPiper.sh"]
