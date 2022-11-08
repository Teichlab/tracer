FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

ARG bowtie_version=2.4.5
ARG samtools_version=1.15.1
ARG trinity_version=v2.14.0
ARG igblast_version=1.19.0
ARG kallisto_version=v0.48.0
ARG salmon_version=1.9.0
ARG gencode_human_version=41
ARG gencode_mouse_version=M30

#Install OS packages
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get -y install wget curl unzip build-essential zlib1g-dev git python3 python3-pip default-jre procps cmake \
    libcairo2-dev pkg-config jellyfish autoconf libbz2-dev liblzma-dev graphviz libgirepository1.0-dev libncurses5-dev \
    python-is-python3 

#Install Bowtie2
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${bowtie_version}/bowtie2-${bowtie_version}-linux-x86_64.zip/download -O /opt/bowtie2-${bowtie_version}-linux-x86_64.zip && \
    cd /opt && \
    unzip bowtie2-${bowtie_version}-linux-x86_64.zip && \
    rm /opt/bowtie2-${bowtie_version}-linux-x86_64.zip

#Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 && \
    tar -xvf samtools-${samtools_version}.tar.bz2 -C /opt && \
    cd /opt/samtools-${samtools_version} && \
    ./configure && \
    make && \
    make install && \
    rm /samtools-${samtools_version}.tar.bz2 

#Install Trinity
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-${trinity_version}/trinityrnaseq-${trinity_version}.FULL.tar.gz && \
    tar -xvzf trinityrnaseq-${trinity_version}.FULL.tar.gz -C /opt && \
    cd /opt/trinityrnaseq-${trinity_version} && \
    make && \
    rm /trinityrnaseq-${trinity_version}.FULL.tar.gz

#Install igblast
RUN wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${igblast_version}/ncbi-igblast-${igblast_version}-x64-linux.tar.gz && \
    tar -xvf ncbi-igblast-${igblast_version}-x64-linux.tar.gz -C /opt && \
    rm /ncbi-igblast-${igblast_version}-x64-linux.tar.gz 

##Set IGDATA variable to see the internal_files
ENV IGDATA="/opt/ncbi-igblast-${igblast_version}"

#Install Kallisto
RUN wget https://github.com/pachterlab/kallisto/releases/download/${kallisto_version}/kallisto_linux-${kallisto_version}.tar.gz && \
    tar -xzvf kallisto_linux-${kallisto_version}.tar.gz -C /opt && \
    rm kallisto_linux-${kallisto_version}.tar.gz

#Install Salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v${salmon_version}/salmon-${salmon_version}_linux_x86_64.tar.gz && \
    tar -xzvf salmon-${salmon_version}_linux_x86_64.tar.gz -C /opt && \
    rm salmon-${salmon_version}_linux_x86_64.tar.gz

#Copy gencode_parse.py from GitHub repo
COPY gencode_parse.py /gencode_parse.py

#Installing transcript sequences for human and mouse 
RUN mkdir -p /var/GRCh38 && \
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${gencode_human_version}/gencode.v${gencode_human_version}.transcripts.fa.gz -O /var/GRCh38/gencode.v${gencode_human_version}.transcripts.fa.gz && \
    gunzip /var/GRCh38/gencode.v${gencode_human_version}.transcripts.fa.gz && \
    python3 /gencode_parse.py /var/GRCh38/gencode.v${gencode_human_version}.transcripts.fa && \
    mv /genemap.tsv /var/GRCh38/ && \
    mv /transcripts.fasta /var/GRCh38/ && \
    rm /var/GRCh38/gencode.v${gencode_human_version}.transcripts.fa

RUN mkdir -p /var/GRCm38 && \
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${gencode_mouse_version}/gencode.v${gencode_mouse_version}.transcripts.fa.gz -O /var/GRCm38/gencode.v${gencode_mouse_version}.transcripts.fa.gz && \
    gunzip /var/GRCm38/gencode.v${gencode_mouse_version}.transcripts.fa.gz && \
    python3 /gencode_parse.py /var/GRCm38/gencode.v${gencode_mouse_version}.transcripts.fa && \
    mv /genemap.tsv /var/GRCm38/ && \
    mv /transcripts.fasta /var/GRCm38/ && \
    rm /var/GRCm38/gencode.v${gencode_mouse_version}.transcripts.fa && \
    rm /gencode_parse.py

#Copy requirments.txt, setup.py and tracerlib from github repo

COPY requirements.txt /requirements.txt

COPY setup.py /setup.py

RUN mkdir -p /opt/tracer

COPY tracerlib /opt/tracer/tracerlib

#Install tracer

RUN pip3 install pycairo

RUN cd /opt/tracer && \ 
    python /setup.py install && \
    rm /setup.py && \
    rm /requirements.txt

#Setting up tracer config file

COPY tracer.conf /home/.tracerrc 

#Adding software to path

ENV PATH="${PATH}:/opt/kallisto:/opt/bowtie2-${bowtie_version}-linux-x86_64:/opt/ncbi-igblast-${igblast_version}/bin:/opt/salmon-${salmon_version}_linux_x86_64/bin:/opt/samtools-${samtools_version}:/opt/trinityrnaseq-${trinity_version}:/opt/tracer"

#Building Indexes
RUN kallisto index -i /var/GRCh38/kallisto.idx /var/GRCh38/transcripts.fasta

RUN kallisto index -i /var/GRCm38/kallisto.idx /var/GRCm38/transcripts.fasta

RUN salmon index --index /var/GRCh38/salmon --transcripts /var/GRCh38/transcripts.fasta

RUN salmon index --index /var/GRCm38/salmon --transcripts /var/GRCm38/transcripts.fasta

#Copying resources

RUN mkdir -p /usr/local/bin/resources

COPY resources /usr/local/bin/resources

#Saving software versions to a file
RUN echo "bowtie2 version: ${bowtie_version}" >> versions.txt && \
    echo "samtools version: ${samtools_version}" >> versions.txt && \
    echo "trinity version: ${trinity_version}" >> versions.txt && \
    echo "igblast version: ${igblast_version}" >> versions.txt && \
    echo "kallisto version: ${kallisto_version}" >> versions.txt && \
    echo "salmon version: ${salmon_version}" >> versions.txt && \
    echo "human transcript version: ${gencode_human_version}" >> versions.txt && \
    echo "mouse transcript version: ${gencode_mouse_version}" >> versions.txt
