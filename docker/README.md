# How to set up Docker

## Don't have Docker? No problem!

Install for **Debian-based** Linux distros:

```bash
sudo apt update && sudo apt install docker.io
```

Install for **CentOS-based** Linux distros:

```bash
## Download rpm repo
sudo yum install -y yum-utils
sudo yum-config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo

## Install Docker engine
sudo yum install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

## Start Docker engine
sudo systemctl start docker
```

Other installations processes [here](https://docs.docker.com/engine/install/)

## Set up the local repository

You need to build your local repository in this way:

```
 .
 |__ envs/
 |__ scripts/
 |__ Dockerfile
 |__ AmpliPiper.sh
```

In this sense, when you run AmpliPiper.sh from within this repository, it will immediately find `envs` and `scripts`.

## Take a look to Dockerfile

What does Dockerfile do?

1. `FROM conda/miniconda3:latest`: Installs conda and python by pulling the latest miniconda image
2. `WORKDIR /app`: Sets the directory where Docker will be looking when running our scripts
3. `ADD . /app`: Adds your local folder to the working directory, which means that it will add AmpliPiper.sh, envs/ and scripts/ (everything that we need!)
4. `ENTRYPOINT ["/bin/bash", "AmpliPiper.sh"]`: Executes AmpliPiper.sh whenever the Docker image is run

## Build the image

You can easily build the image from the Dockerfile: you just have to move to the directory where it is placed and run:

```bash
docker build . -t nhmvienna/amplipiper:latest
```

Let's quickly break down the command:

- `build` looks for a Dockerfile in the '.' directory, and executes all the commands inside it
- `-t` flag stands for 'tag': as you can see, the tag has three portions:
    
    + the owner of the image (I put nhmvienna as an example, **but this portion is the only one that has to coincide with a registered Docker account**)
    + the name of the image (amplipiper)
    + the actual tag (or version, which is 'latest')

## Push the image to Docker hub

Once the build is complete, you may want to push the image to Docker hub:

```bash
docker push nhmvienna/amplipiper:latest
```

Depending on your internet connection and on the size of the image, this may take a while!

## Run the image

If you were a new user and you did not build the image from scratch, you now would have to first pull the image from the hub:

```bash
docker pull nhmvienna/amplipiper:latest
```

And the run it:

```bash
docker run nhmvienna/amplipiper:latest [OPTIONS] 
```

You can pass all the options of the pipeline directly to the Docker run command!

If you wish to make the pipeline run in the background, you can also add the `-d` flag:

```bash
docker run -d nhmvienna/amplipiper:latest [OPTIONS] 
```

And you're done!
