FROM conda/miniconda3:latest 

## Your local directory should be structured like this:
## .
## |__ envs/
## |__ scripts/
## |__ Dockerfile
## |__ AmpliPiper.sh

WORKDIR /app

ADD . /app

ENTRYPOINT ["/bin/bash", "AmpliPiper.sh"]
