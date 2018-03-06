Snakemake molecular trait pipeline
==================================

Snakemake pipeline for handling molecular traits.

* GitHub repo: https://github.com/letaylor/snakemake-qtl
* Free software: MIT license


Overview
--------

The workflow is divided into modules. Each module contains a Snakemake file that will execute that analysis chunk. All code specific to a module is housed within the module's directory (e.g., scripts and wrappers). General code that may be used by multiple modules is housed in this directory. In addition, this directory houses the data subdirectory which contains data for minimal examples that may be referenced by multiple modules. Each module is tracked with a version at the top of the Snakemake file and in the setup.cfg file.

1. Latent factor analysis using PEER

2. QTL mapping


Quick Setup
-----------

A quick demo is provided for each module. For instance, see the README file in the qtl module. 

```bash
# clone the repo
git clone https://github.com/letaylor/snakemake-qtl

# set code base path
SNK_REPO="`pwd`/snakemake-qtl"

# run README demo in the qtl module
# alternatively keep reading for reproducible setup
```

If you are running jobs on the cluster, it is best to first start a [tmux](https://github.com/tmux/tmux) session so that the session may be re-attached at a later time point. 

```bash
# start session
tmux new -s snkmk

# log back into a session
tmux -a snkmk
```


Reproducibility: Conda   
----------------------

One option to enhance reproducibility is to install software used via Conda. Note that some software (e.g., QTLtools) will still be missing since they are not currently in a Conda channel.

```bash
# install environment using Conda
conda env create --name moltraits --file ${SNK_REPO}/environment.yml

# activate the new Conda environment
source activate moltraits
```


Reproducibility: Docker
-----------------------

An alternative option is to use a Docker image. One can easily generate a Docker image of all of the software used in this repository (including software missing from a Conda channel) by using the Dockerfile. 

```bash
# build the Docker image using the Dockerfile
cd ${SNK_REPO}
docker build -t moltraits .


# build the Docker image using jupyter-repo2docker 
# here, one would delete the current Dockerfile 
# and add the additional install files via postBuild
jupyter-repo2docker --no-build --debug https://github.com/letaylor/snakemake-qtl
```

A pre-compiled Docker image is housed on the Docker Cloud. Download and use this Docker image:

```bash
# download the Docker image 
docker pull letaylor/moltraits:latest

# run the Snakemake pipeline through the container
# the Snakemake command used below is described in the qtl/README.md file
docker run -t -v ${SNK_REPO}:/SNK_REPO -v $(pwd):/CUR_DIR -e USERID=$UID letaylor/moltraits:latest "snakemake --snakefile /SNK_REPO/qtl/Snakefile --directory /CUR_DIR --configfile /CUR_DIR/config_analysis.json --printshellcmds"
```


Reproducibility: Singularity
----------------------------

A final option is to load the above Docker image using Singularity, which is designed for high-performance compute systems. To do so, simply add the --use-singularity flag when calling snakemake as descibed in the other README.md docs (within the different modules).

As an example, see below. Note that the Docker image is specified in the DOCKER variable of the config file (config_analysis.json).

```bash
# the Snakemake command used below is described in the qtl/README.md file
snakemake --snakefile ${SNK_REPO}/qtl/Snakefile --configfile config_analysis.json --printshellcmds --use-singularity --singularity-prefix $(pwd)/.snakemake/singularity
```
