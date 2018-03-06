FROM continuumio/miniconda3:latest

# use bash as default shell, rather than sh
ENV SHELL /bin/bash

# set up a user, to avoid running containers as root
ENV NB_USER container_user
ENV HOME /home/${NB_USER}
RUN adduser --disabled-password \
    --gecos "Default user" \
    ${NB_USER}

# download things to tmp
WORKDIR /tmp

# download and install the conda environment
ENV ENV_FILE https://raw.githubusercontent.com/letaylor/snakemake-qtl/master/environment.yml
RUN wget ${ENV_FILE} \
    && conda update conda --yes \
    && conda env create --name container_env --file 'environment.yml' \
    && conda clean --all --yes \
    && conda list --name container_env \
    #&& echo 'source activate container_env' >> .bashrc \
    && rm 'environment.yml'

# export the conda path, since .bashrc is not automatically sourced 
# when the image is used
ENV PATH /opt/conda/envs/container_env/bin:$PATH
ENV CONDA_DEFAULT_ENV container_env
ENV CONDA_PREFIX /opt/conda/envs/container_env

# install QTLTools
RUN wget \
https://qtltools.github.io/qtltools/binaries/QTLtools_1.1_Ubuntu14.04_x86_64.tar.gz \
    && tar -zxvf \
        QTLtools_1.1_Ubuntu14.04_x86_64.tar.gz \
        QTLtools_1.1_Ubuntu14.04_x86_64 \
    && rm QTLtools_1.1_Ubuntu14.04_x86_64.tar.gz \
    && chmod 755 QTLtools_1.1_Ubuntu14.04_x86_64 \
    && mv QTLtools_1.1_Ubuntu14.04_x86_64 /bin/QTLtools

# install veqtl-mapper
RUN wget https://github.com/funpopgen/veqtl-mapper/releases/download/1.02/veqtl-mapper \
    && chmod 755 veqtl-mapper \
    && mv veqtl-mapper /bin/veqtl-mapper

# set wd to user home
WORKDIR ${HOME}

# clear tmp if there is anything in there... 
RUN rm -rf /tmp/*

# set the user
USER ${NB_USER}
