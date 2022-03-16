FROM continuumio/miniconda3:latest

LABEL maintainer="matt.demaere@gmail.com"
LABEL org.label-schema.schema-version="1.0"
LABEL org.label-schema.name="cerebis/bin3c"
LABEL org.label-schema.description="bin3C - extract metagenome-assembled genomes (MAGs) from metagenomic data using Hi-C"
LABEL org.label-schema.url="http://github.com/cerebis/bin3C/"
LABEL org.label-schema.vcs-url="http://github.com/cerebis/bin3C/"
LABEL org.label-schema.vcs-ref="7425ea62b9d90ddee405b46baa2806e6290ca3ce"
LABEL org.label-schema.version="0.4"
LABEL org.label-schema.docker.cmd="docker run -v /path/to/data:/app cerebis/bin3c --help"

RUN apt-get update && \
    apt-get install -y apt-utils && \
    apt-get install -y \
        make \
        gcc \
        g++ \
        git \
        llvm \
        libncurses-dev \
        libblas-dev \
        wget \
        bzip2 \
        gzip \
        xz-utils \
        zlib1g-dev && \
    apt-get clean

RUN conda install --yes -c conda-forge -c bioconda "python==3.7" bwa samtools && \
    conda clean -afy

RUN pip3 install --no-cache-dir cython numpy && \
    pip3 install --no-cache-dir git+https://github.com/cerebis/bin3C@py3

RUN mkdir -p /app
WORKDIR /app
ENTRYPOINT ["bin3C"]
CMD ["--help"]

