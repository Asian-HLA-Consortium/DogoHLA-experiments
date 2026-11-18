FROM ubuntu:24.04

LABEL org.opencontainers.image.authors="Asian-HLA-Consortium / DogoHLA Team"
LABEL org.opencontainers.image.description="DogoHLA typing and panel workflow container"
LABEL org.opencontainers.image.version="1.0.0"

ENV LC_ALL=C
ENV LANG=C
ENV PATH="/opt/t1k:/opt/minimap2:/opt/samtools/bin:/opt/jellyfish/bin:/opt/bwa:${PATH}"
ENV PYTHONPATH="/opt/immuannot:${PYTHONPATH}"

# System dependencies and build toolchain
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        curl \
        wget \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libncurses5-dev \
        python3 \
        python3-dev \
        python3-pip \
        python3-setuptools \
        bwa \
        samtools \
        minimap2 \
        jellyfish \
        ca-certificates && \
    pip install --no-cache-dir --break-system-packages \
        "numpy==2.5.4" \
        "scipy==1.18.1" \
        "pandas==3.1.0" \
        "pysam==0.24.1" \
        "edlib==1.3.9.post1" && \
    mkdir -p /opt/t1k && cd /opt/t1k && \
    wget -q https://github.com/mourisl/T1K/archive/refs/tags/v1.0.6.tar.gz -O t1k.tar.gz && \
    tar -xzf t1k.tar.gz --strip-components=1 && \
    make -j4 && \
    rm -f t1k.tar.gz && \
    mkdir -p /opt/immuannot && cd /opt/immuannot && \
    git clone https://github.com/YingZhou001/Immuannot.git . && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["/bin/bash"]
