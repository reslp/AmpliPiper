# Technical documentation for AmpliPiper Docker image

The right context to build AmpliPiper Docker image can be found [here](https://github.com/AstraBert/amplipiper)

## Dockerfile

Based on [`conda/miniconda3`](https://hub.docker.com/r/conda/miniconda3/) and on [`mambaorg/micromamba`](https://hub.docker.com/r/mambaorg/micromamba) Docker images, the Dockerfile includes a:

- Multi-staged build: `mambaorg/micromamba` is mounted as `mambabase` and `micromamba` is copied from there into `/usr/bin` of our base image (`conda/miniconda3`)
- Working directory is set as `/app/`
- Context transferring from local filesystem to Docker image
- Execution of the setup script, which has been modified to match with `micromamba` availability instead of `mamba` (find it at [`setup.docker.sh`](./setup.docker.sh))

## CI/CD

Continuous-Integration/Continous-Developement is maintained through [the base repo for this Docker image](https://github.com/AstraBert/amplipiper) and managed with [**GitHub Container Registry**](https://ghcr.io/). the image will be availble under: `ghcr.io/astrabert/amplipiper` and we advise you pull the `main` tag, as it is always on track with the modifications operated on the `main` branch of the repository (the one on which the image is based).

## Build locally

Execute the code to build the image locally:

```bash
git clone https://github.com/AstraBert/amplipiper.git
cd amplipiper
docker build . -t YOUR-USER-NAME/IMAGE-NAME:tag
```