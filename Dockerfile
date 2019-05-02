# run as:
#   docker build -t tracer .

#start off with a plain Debian
FROM debian:latest

#basic setup stuff, including bowtie2
RUN apt-get update && apt-get -y upgrade
RUN apt-get -y install wget curl unzip build-essential zlib1g-dev git python3 python3-pip bowtie2 openjdk-8-jre procps

#Trinity - depends on zlib1g-dev and openjdk-8-jre installed previously
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.4.0.zip
RUN unzip Trinity-v2.4.0.zip && rm Trinity-v2.4.0.zip
RUN cd /trinityrnaseq-Trinity-v2.4.0 && make

#IgBLAST, plus the setup of its super weird internal_data thing. don't ask. just needs to happen
#and then on top of that, the environmental variable thing facilitates the creation of a shell wrapper. fun
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.7.0/ncbi-igblast-1.7.0-x64-linux.tar.gz
RUN tar -xzvf ncbi-igblast-1.7.0-x64-linux.tar.gz && rm ncbi-igblast-1.7.0-x64-linux.tar.gz
RUN cd /ncbi-igblast-1.7.0/bin/ && wget -r ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data && \
	wget -r ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/optional_file && \
	mv ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data . && \
	mv ftp.ncbi.nih.gov/blast/executables/igblast/release/optional_file . && \
	rm -r ftp.ncbi.nih.gov

#aligners - kallisto and salmon
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz
RUN tar -xzvf kallisto_linux-v0.43.1.tar.gz && rm kallisto_linux-v0.43.1.tar.gz
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.8.2/Salmon-0.8.2_linux_x86_64.tar.gz
RUN tar -xzvf Salmon-0.8.2_linux_x86_64.tar.gz && rm Salmon-0.8.2_linux_x86_64.tar.gz

#graphviz, which lives in a sufficient form (dot/neato) in apt-get apparently
RUN apt-get -y install graphviz

#tracer proper
COPY . /tracer
RUN cd /tracer && pip3 install -r docker_helper_files/requirements_stable.txt && python3 setup.py install

#obtaining the transcript sequences. no salmon/kallisto indices as they make dockerhub unhappy for some reason
RUN mkdir GRCh38 && cd GRCh38 && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz && \
	gunzip gencode.v27.transcripts.fa.gz && python3 /tracer/docker_helper_files/gencode_parse.py gencode.v27.transcripts.fa && rm gencode.v27.transcripts.fa
RUN mkdir GRCm38 && cd GRCm38 && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M15/gencode.vM15.transcripts.fa.gz && \
	gunzip gencode.vM15.transcripts.fa.gz && python3 /tracer/docker_helper_files/gencode_parse.py gencode.vM15.transcripts.fa && rm gencode.vM15.transcripts.fa

#placing a preconfigured tracer.conf in ~/.tracerrc
RUN cp /tracer/docker_helper_files/docker_tracer.conf ~/.tracerrc

#this is a tracer container, so let's point it at a tracer wrapper that sets the silly IgBLAST environment variable thing
ENTRYPOINT ["bash", "/tracer/docker_helper_files/docker_wrapper.sh"]
