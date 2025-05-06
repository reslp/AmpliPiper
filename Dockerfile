FROM reslp/mamba:2.0.5

WORKDIR /app

COPY . /app/AmpliPiper/
#RUN git clone https://github.com/reslp/AmpliPiper.git

RUN bash /app/AmpliPiper/shell/setup.sh

RUN chmod +x /app/AmpliPiper/shell/AmpliPiper.sh

ENV PATH=$PATH:/app/AmpliPiper/shell
#RUN sed -i 's|tmp=$(dirname $0)|tmp=/app/AmpliPiper/shell|' /app/AmpliPiper/shell/AmpliPiper.sh
#RUN apt update && apt install -y vim

ENV AMICONTAINER="yes"


CMD ["AmpliPiper.sh"]
