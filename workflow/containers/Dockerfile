FROM ubuntu:24.04

LABEL org.opencontainers.image.authors="Asian-HLA-Consortium / DogoHLA Team"
LABEL org.opencontainers.image.description="DogoHLA typing and panel workflow container"
LABEL org.opencontainers.image.version="1.0.0"

ENV LC_ALL=C
ENV LANG=C
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/immuannot/scripts.pub.v3:/opt/t1k:/usr/local/bin:${PATH}"
ENV PYTHONPATH="/opt/immuannot/scripts.pub.v3:/opt/immuannot:${PYTHONPATH}"
ENV IMMUANNOT_SH="/opt/immuannot/scripts.pub.v3/immuannot.sh"
ENV IMMUANNOT_REF="/opt/immuannot/Data-ipd3590-kir2130"

# 1. System dependencies and build toolchain
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        curl \
        wget \
        unzip \
        pkg-config \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libncurses-dev \
        python3 \
        python3-dev \
        python3-pip \
        python3-setuptools \
        ca-certificates && \
    # 2. Python packages (exact pinned versions)
    pip install --no-cache-dir --break-system-packages \
        "numpy==2.5.4" \
        "scipy==1.18.1" \
        "pandas==3.1.0" \
        "pysam==0.24.1" \
        "edlib==1.3.9.post1" && \
    # 3. Source builds for tools differing from Ubuntu 24.04 apt versions:
    # 3a. bwa 0.7.19 (0.7.19-r1273)
    mkdir -p /tmp/build-bwa && cd /tmp/build-bwa && \
    wget -q https://github.com/lh3/bwa/archive/refs/tags/v0.7.19.tar.gz -O bwa.tar.gz && \
    tar -xzf bwa.tar.gz --strip-components=1 && \
    make -j4 && \
    cp bwa /usr/local/bin/ && \
    rm -rf /tmp/build-bwa && \
    # 3b. samtools 1.24 (HTSlib 1.24; apt has 1.19.2)
    mkdir -p /tmp/build-samtools && cd /tmp/build-samtools && \
    wget -q https://github.com/samtools/samtools/releases/download/1.24/samtools-1.24.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjf samtools.tar.bz2 --strip-components=1 && \
    ./configure --prefix=/usr/local && \
    make -j4 && \
    make install && \
    rm -rf /tmp/build-samtools && \
    # 3c. minimap2 2.31 (2.31-r1302; apt has 2.26)
    mkdir -p /tmp/build-minimap2 && cd /tmp/build-minimap2 && \
    wget -q https://github.com/lh3/minimap2/releases/download/v2.31/minimap2-2.31.tar.bz2 -O minimap2.tar.bz2 && \
    tar -xjf minimap2.tar.bz2 --strip-components=1 && \
    make -j4 && \
    cp minimap2 /usr/local/bin/ && \
    rm -rf /tmp/build-minimap2 && \
    # 3d. jellyfish 2.2.10 (apt has 2.3.1)
    mkdir -p /tmp/build-jellyfish && cd /tmp/build-jellyfish && \
    wget -q https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz -O jellyfish.tar.gz && \
    tar -xzf jellyfish.tar.gz --strip-components=1 && \
    ./configure --prefix=/usr/local && \
    make -j4 && \
    make install && \
    ldconfig && \
    rm -rf /tmp/build-jellyfish && \
    # 4. T1K v1.0.6 (v1.0.6-r216)
    mkdir -p /opt/t1k && cd /opt/t1k && \
    wget -q https://github.com/mourisl/T1K/archive/refs/tags/v1.0.6.tar.gz -O t1k.tar.gz && \
    tar -xzf t1k.tar.gz --strip-components=1 && \
    make -j4 && \
    ln -s /opt/t1k/run-t1k /usr/local/bin/run-t1k && \
    ln -s /opt/t1k/fastq-extractor /usr/local/bin/fastq-extractor && \
    rm -f t1k.tar.gz && \
    # 5. Immuannot (commit 31362c370123f6706a04ffdda6a3bf0c108328e0) and reference bundle Data-ipd3590-kir2130
    mkdir -p /opt/immuannot && cd /opt/immuannot && \
    git clone https://github.com/YingZhou001/Immuannot.git . && \
    git checkout 31362c370123f6706a04ffdda6a3bf0c108328e0 && \
    mkdir -p /opt/immuannot/ipd && cd /opt/immuannot/ipd && \
    curl -sfL -o hla.dat.zip https://raw.githubusercontent.com/ANHIG/IMGTHLA/3590/hla.dat.zip && \
    unzip -o -q hla.dat.zip && \
    mv hla.dat hla.3.59.0.dat && \
    curl -sfL -o kir.2.13.0.dat https://raw.githubusercontent.com/ANHIG/IPDKIR/2130/kir.dat && \
    P=/opt/immuannot/scripts.pub.v3/prepare-reference && \
    REF=/opt/immuannot/Data-ipd3590-kir2130 && \
    mkdir -p "$REF" && \
    mkdir -p /tmp/predata && cd /tmp/predata && \
    bash /opt/immuannot/scripts.pub.v3/predata.sh /opt/immuannot/ipd/hla.3.59.0.dat /opt/immuannot/ipd/kir.2.13.0.dat "$P/c4a.gb" "$REF" && \
    cat "$P/curated/DRB3.ele.csv" | gzip >> "$REF/alleles.csv.gz" && \
    cat "$P/curated/DRB3.new.fa" | gzip >> "$REF/gen.fa.gz" && \
    cat "$P/curated/DRB3.new.cds.fa" | gzip >> "$REF/CDSseq/HLA-DRB3.fa.gz" && \
    rm -rf /tmp/predata /opt/immuannot/ipd/hla.dat.zip && \
    ln -s /opt/immuannot/scripts.pub.v3/immuannot.sh /usr/local/bin/immuannot.sh && \
    # Cleanup apt caches
    apt-get clean && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["/bin/bash"]
