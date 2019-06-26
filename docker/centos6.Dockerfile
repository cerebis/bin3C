#
# stage 1 - build
#
FROM centos:centos6 as builder
MAINTAINER Matthew DeMaere "matt.demaere@gmail.com"

# create up to date system
RUN yum update -y && \
    yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-6.noarch.rpm && \
    yum install -y \
        centos-release-scl \
        bzip2 \
        bzip2-libs \
        bzip2-devel \
        freetype-devel \
        git \
        libcurl-devel \
        libpng-devel \
        libgomp \
        llvm \
        ncurses-devel \
        openblas-devel \
        openssl-devel \
        python2 \
        python2-devel \
        python2-pip \
        redhat-rpm-config \
        wget \
        xz-devel \
        zlib-devel && \
    yum install -y python27 devtoolset-6-gcc devtoolset-6-gcc-c++ devtoolset-6-gcc-gfortran && \
    yum clean all -y

# setup more current python and gcc support
SHELL ["/usr/bin/scl", "enable",  "python27", "devtoolset-6"]

# first update pip and install cython
RUN pip2 install -U pip setuptools && pip2 install --user cython

# install qc3C pgtk branch
RUN pip2 install --user "numpy<1.15" && \
    pip2 install --user git+https://github.com/cerebis/bin3C@pgtk

# install spades
WORKDIR /opt/
RUN wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz && \
    tar xzf SPAdes-3.13.0-Linux.tar.gz

RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar xjf bwa-0.7.17.tar.bz2
WORKDIR bwa-0.7.17
RUN make

WORKDIR /usr/local
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar xjf htslib-1.9.tar.bz2
WORKDIR /usr/local/htslib-1.9
RUN ./configure && make

WORKDIR /usr/local/
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xjf samtools-1.9.tar.bz2
WORKDIR /usr/local/samtools-1.9
RUN ./configure && make

WORKDIR /opt
RUN wget -O bbmap.tar.gz "https://downloads.sourceforge.net/project/bbmap/BBMap_38.44.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2Ffiles%2Flatest%2Fdownload&ts=1555391525" && \
    tar xzf bbmap.tar.gz

#
# stage two - runtime
#
FROM centos:centos6

# restore minimum runtime dependencies
RUN yum update -y && \
    yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-6.noarch.rpm && \
    yum install -y \
        centos-release-scl \
        bzip2 \
        freetype \
        java-1.8.0-openjdk \
        libpng \
        libgomp \
        llvm \
        openblas \
        python2 \
        xz && \
    yum install -y python27 && \
    yum clean all -y && \
    ldconfig

# copy over the installed python packages
COPY --from=builder /opt/rh/python27/root /opt/rh/python27/root/
COPY --from=builder /root/.local/bin/bin3C /opt/rh/python27/root/usr/bin/
COPY --from=builder /root/.local/lib/python2.7/site-packages /opt/rh/python27/root/usr/lib/python2.7/site-packages/
COPY --from=builder /usr/local/samtools-1.9/samtools /usr/local/bin/
COPY --from=builder /usr/local/htslib-1.9/bgzip /usr/local/bin/
COPY --from=builder /opt/bwa-0.7.17/bwa /usr/local/bin/
COPY --from=builder /opt/bbmap /opt/bbmap/
COPY --from=builder /opt/SPAdes-3.13.0-Linux /opt/SPAdes-3.13.0-Linux/
RUN ln -s /opt/bbmap/*sh /usr/local/bin/ && \
    ln -s /opt/SPAdes-3.13.0-Linux/bin/* /usr/local/bin/ && \
    chmod 755 /opt/rh/python27/root/usr/lib/python2.7/site-packages

# further setup
COPY ./root/ /

ENV HOME=/opt/app-root/ \
    PATH=/opt/app-root/bin:/opt/rh/python27/root/usr/bin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN mkdir -p ${HOME}

WORKDIR ${HOME}

ENV BASH_ENV=/opt/app-root/etc/scl_enable \
    ENV=/opt/app-root/etc/scl_enable \
    PROMPT_COMMAND=". /opt/app-root/etc/scl_enable"

# Set the default CMD to print the usage of the language image
ENTRYPOINT ["/usr/bin/container-entrypoint"]
CMD ["usage"]
