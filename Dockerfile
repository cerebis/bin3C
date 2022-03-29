#
# stage one - build
#
FROM continuumio/miniconda3:latest AS builder

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

# create base conda environment
COPY environment.yaml .
RUN conda env create -f environment.yaml 
# ensure RUN commands happen within an activated environment
RUN echo "conda activate bin3c_env" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]
# pip install bin3C
RUN pip install --no-cache-dir git+https://github.com/cerebis/bin3C@py3

# pack up the environment and make it a stand-alone directory
RUN conda install -c conda-forge conda-pack
RUN conda-pack -n bin3c_env -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && rm /tmp/env.tar
RUN /venv/bin/conda-unpack
 
#
# stage two - runtime
#
FROM debian:buster AS runtime

LABEL maintainer="matt.demaere@gmail.com"
LABEL org.label-schema.schema-version="1.0"
LABEL org.label-schema.name="cerebis/bin3c"
LABEL org.label-schema.description="bin3C - extract metagenome-assembled genomes (MAGs) from metagenomic data using Hi-C"
LABEL org.label-schema.url="http://github.com/cerebis/bin3C/"
LABEL org.label-schema.vcs-url="http://github.com/cerebis/bin3C/"
LABEL org.label-schema.vcs-ref="7425ea62b9d90ddee405b46baa2806e6290ca3ce"
LABEL org.label-schema.version="0.4"
LABEL org.label-schema.docker.cmd="docker run -v /path/to/data:/app cerebis/bin3c"

# copy the stand-alone environment to runtime image
COPY --from=builder /venv /venv
COPY entrypoint.sh /venv/bin/

ENTRYPOINT ["/venv/bin/entrypoint.sh", "bin3C"]
CMD ["--help"]

