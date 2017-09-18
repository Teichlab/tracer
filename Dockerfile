# run as:
#   docker build -t tracer .

#start off with a plain Debian
FROM debian:latest

#basic setup stuff, including bowtie2
RUN apt-get update && apt-get -y upgrade
RUN apt-get -y install wget curl unzip build-essential zlib1g-dev git python3 python3-pip bowtie2 openjdk-8-jre

#Trinity - depends on zlib1g-dev and openjdk-8-jre installed previously
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.4.0.zip
RUN unzip Trinity-v2.4.0.zip && rm Trinity-v2.4.0.zip
RUN cd /trinityrnaseq-Trinity-v2.4.0 && make

#IgBLAST, plus the setup of its super weird internal_data thing. don't ask. just needs to happen
#and then on top of that, the environmental variable thing facilitates the creation of a shell wrapper. fun
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.7.0/ncbi-igblast-1.7.0-x64-linux.tar.gz
RUN tar -xzvf ncbi-igblast-1.7.0-x64-linux.tar.gz && rm ncbi-igblast-1.7.0-x64-linux.tar.gz
RUN cd /ncbi-igblast-1.7.0/bin/ && wget -r ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data && \
	mv ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data . && rm -r ftp.ncbi.nih.gov

#aligners - kallisto and salmon
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz
RUN tar -xzvf kallisto_linux-v0.43.1.tar.gz && rm kallisto_linux-v0.43.1.tar.gz
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.8.2/Salmon-0.8.2_linux_x86_64.tar.gz
RUN tar -xzvf Salmon-0.8.2_linux_x86_64.tar.gz && rm Salmon-0.8.2_linux_x86_64.tar.gz

#graphviz, along with its sea of dependencies that otherwise trip up the dpkg -i
RUN apt-get -y install libgd3 libgts-0.7-5 liblasi0 libltdl7 freeglut3 libglade2-0 libglu1-mesa libglu1 libgtkglext1 libxaw7
RUN wget http://www.graphviz.org/pub/graphviz/stable/ubuntu/ub13.10/x86_64/libgraphviz4_2.38.0-1~saucy_amd64.deb
RUN dpkg -i libgraphviz4_2.38.0-1~saucy_amd64.deb && apt-get -y -f install
RUN wget http://www.graphviz.org/pub/graphviz/stable/ubuntu/ub13.10/x86_64/graphviz_2.38.0-1~saucy_amd64.deb
RUN dpkg -i graphviz_2.38.0-1~saucy_amd64.deb && apt-get -y -f install
RUN rm libgraphviz4_2.38.0-1~saucy_amd64.deb && rm graphviz_2.38.0-1~saucy_amd64.deb

#tracer proper, along with repositioning its test_data and resources in a manner the installed package can understand
COPY . /tracer
RUN cd /tracer && pip3 install -r requirements.txt && python3 setup.py install
RUN cp -r /tracer/test_data /usr/local/lib/python3.5/dist-packages/tracer-0.5-py3.5.egg/
RUN cp -r /tracer/resources /usr/local/lib/python3.5/dist-packages/tracer-0.5-py3.5.egg/

#obtaining the transcript sequences and kallisto/salmon indices
RUN mkdir GRCh38 && cd GRCh38 && wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz && \
	gunzip gencode.v27.transcripts.fa.gz && python3 /tracer/docker_helper_files/gencode_parse.py gencode.v27.transcripts.fa && rm gencode.v27.transcripts.fa && \
	/Salmon-0.8.2_linux_x86_64/bin/salmon index -t transcripts.fasta -i salmon && \
	/kallisto_linux-v0.43.1/kallisto index -i kallisto.idx transcripts.fasta
RUN mkdir GRCm38 && cd GRCm38 && wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M15/gencode.vM15.transcripts.fa.gz && \
	gunzip gencode.vM15.transcripts.fa.gz && python3 /tracer/docker_helper_files/gencode_parse.py gencode.vM15.transcripts.fa && rm gencode.vM15.transcripts.fa && \
	/Salmon-0.8.2_linux_x86_64/bin/salmon index -t transcripts.fasta -i salmon && \
	/kallisto_linux-v0.43.1/kallisto index -i kallisto.idx transcripts.fasta

#world's longest sed pipe, turning the github demo tracer.conf into a functional file representing where things live
RUN sed 's/#igblastn_path = \/path\/to\/igblastn/igblastn_path = \/ncbi-igblast-1.7.0\/bin\/igblastn/g' /tracer/tracer.conf | \
	sed 's/#makeblastdb_path = \/path\/to\/makeblastdb/makeblastdb_path = \/ncbi-igblast-1.7.0\/bin\/makeblastdb/g' | \
	sed 's/#kallisto_path = \/path\/to\/kallisto/kallisto_path = \/kallisto_linux-v0.43.1\/kallisto/g' | \
	sed 's/#salmon_path = \/path\/to\/salmon/salmon_path = \/Salmon-0.8.2_linux_x86_64\/bin\/salmon/g' | \
	sed 's/#trinity_path = \/path\/to\/trinity/trinity_path = \/trinityrnaseq-Trinity-v2.4.0\/Trinity/g' | \
	sed ':a;N;$!ba;s/\[base_transcriptomes\]\n# reference transcriptomes for kallisto\/salmon\nMmus = \/path\/to\/kallisto\/transcriptome_for_Mmus\nHsap = \/path\/to\/kallisto\/transcriptome_for_Hsap/\[base_transcriptomes\]\n# reference transcriptomes for kallisto\/salmon\nMmus = \/GRCm38\/transcripts.fasta\nHsap = \/GRCh38\/transcripts.fasta/g' | \
	sed ':a;N;$!ba;s/\[salmon_base_indices\]\n# salmon indices created from \[base_transcriptomes\] above; needed only when option --small_index is used\nMmus = \/path\/to\/salmon\/index_for_Mmus\nHsap = \/path\/to\/salmon\/index_for_Hsap/\[salmon_base_indices\]\n# salmon indices created from \[base_transcriptomes\] above; needed only when option --small_index is used\nMmus = \/GRCm38\/salmon\nHsap = \/GRCh38\/salmon/g' | \
	sed ':a;N;$!ba;s/\[kallisto_base_indices\]\n# kallisto indices created from \[base_transcriptomes\] above; needed only when option --small_index is used\nMmus = \/path\/to\/kallisto\/index_for_Mmus\nHsap = \/path\/to\/kallisto\/index_for_Hsap/\[kallisto_base_indices\]\n# kallisto indices created from \[base_transcriptomes\] above; needed only when option --small_index is used\nMmus = \/GRCm38\/kallisto.idx\nHsap = \/GRCh38\/kallisto.idx/g' > ~/.tracerrc

#this is a tracer container, so let's point it at a tracer wrapper that sets the silly IgBLAST environment variable thing
ENTRYPOINT ["bash", "/tracer/docker_helper_files/docker_wrapper.sh"]
