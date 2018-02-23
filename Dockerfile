FROM continuumio/miniconda3

# set bash as the default entry point
ENTRYPOINT [ '/bin/bash', '-c' ]

WORKDIR /src

# download and install the conda environment
RUN wget \
https://raw.githubusercontent.com/letaylor/snakemake-qtl/master/environment.yml \
    && conda env create --name moltraits --file environment.yml \
    && echo 'source activate moltraits' > /root/.bashrc

# activate conda moltraits environment
#RUN /bin/bash -c 'source activate moltraits'

# export the conda path, since .bashrc is not automatically sourced 
# when the image is used
ENV PATH /opt/conda/envs/moltraits/bin:$PATH
ENV CONDA_DEFAULT_ENV moltraits
ENV CONDA_PREFIX /opt/conda/envs/moltraits

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
RUN wget https://github.com/funpopgen/veqtl-mapper/releases/download/1.01/veqtl-mapper \
    && chmod 755 veqtl-mapper \
    && mv veqtl-mapper /bin/veqtl-mapper


WORKDIR /



