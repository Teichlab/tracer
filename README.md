# TraCeR
TraCeR - reconstruction of T cell receptor sequences from single-cell RNA-seq data.

##Introduction
This tool reconstructs the sequences of rearranged and expressed T cell receptor genes from single-cell RNA-seq data. It then uses the TCR sequences to identify cells that have the same receptor sequences and so derive from the same original clonally-expanded cell. 

For more information on TraCeR, its validation and how it can be applied to investigate T cell populations during infection, see our publication at *link to MS/publication here*.

##Installation
TraCeR is written in Python and so can just be downloaded, made executable (with `chmod u+x tracer`) and run. Download the latest version and accompanying files from www.github.com/teichlab/tracer. 

Tracer relies on several additional tools and Python modules that you should install.

###Pre-requisites

####Software####
1. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - required for alignment of reads to synthetic TCR genomes.
2. [Trinity](http://trinityrnaseq.github.io) - required for assembly of reads into TCR contigs.
3. [IgBLAST](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone) ([FTP site](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/)) - required for analysis of assembled contigs.
4. [Kallisto](http://pachterlab.github.io/kallisto/) - required for quantification of TCR expression.
5. [Graphviz](http://www.graphviz.org) - Dot and Neato drawing programs required for visualisation of clonotype graphs.

####Python modules####
1. [Matplotlib](http://matplotlib.org)
2. [Seaborn](http://stanford.edu/~mwaskom/software/seaborn/)
3. [Biopython](http://biopython.org/)
4. [Prettytable](https://code.google.com/p/prettytable/)
5. [Levenshtein](https://pypi.python.org/pypi/python-Levenshtein/)
6. [Networkx](https://networkx.github.io)


##Setup##
Once the prerequisites above are installed and working you're ready to tell TraCeR where to find them.

TraCeR uses a configuration file to point it to the locations of files that it needs and a couple of other options. By default, this is `tracer.conf` in the same directory as the TraCeR executable. The `-c` option to the various tracer modules allows you to specify any other file to act as the configuration file. 

###External tool locations###
Edit `tracer.conf` (or a copy) so that the paths within the `[tool_locations]` section point to the executables for all of the required tools. 

		[tool_locations]
		#paths to tools used by TraCeR for alignment, quantitation, etc
		bowtie2_path = /path/to/bowtie2
		igblast_path = /path/to/igblastn
		kallisto_path = /path/to/kallisto
		trinity_path = /path/to/trinity
		dot_path = /path/to/dot
		neato_path = /path/to/neato
		
###Resource locations and necessary files###
The tools used by TraCeR need a variety of additional files to work properly and to allow extraction of TCR-derived reads and expression quantification etc. The locations of these files are specified in the other sections of the configuration file and are detailed below.

Currently, organism-specific files (TCR gene sequences, synthetic genome indices, igblast_indices) for mouse and human are distributed with the source-code in the `resources` directory. There will soon be a `build` module that constructs all of the necessary synthetic genomes and indices from any collection of V and J gene sequences. 

####Bowtie synthetic genomes path####
		[bowtie2_options]
		synthetic_genome_index_path = resources/synthetic_genomes/mouse

This path specifies the directory that contains Bowtie2 indices constructed from all possible combinations of V and J segments for each locus. 

####Trinity options####
#####Jellyfish memory#####
		[trinity_options]
		#line below specifies maximum memory for Trinity Jellyfish component. Set it appropriately for your environment.
		max_jellyfish_memory = 1G

Trinity needs to know the maximum memory available to it for the Jellyfish component. Specify this here. 
#####HPC configuration#####
    trinity_grid_conf = /nfs/research2/teichmann/mike/TCR/scripts/TCR_Trinity.conf

Trinity can parallelise contig assembly by submitting jobs across a compute cluster. If you're running in such an environment you can specify an optional trinity config file here. See the Trinity documentation for more information.

####IgBLAST options####
#####Databases path#####
		[IgBlast_options]
		igblast_index_location = resources/igblast_dbs/mouse
#####VDJ sequences#####
This path specifies the directory that contains IgBLAST database files for V, D and J genes. These files are named `imgt_tcr_db_<SEGMENT>.fa`.

    imgt_seq_location = resources/imgt_sequences/mouse
		
Path to fasta files with sequences for each V, D or J gene. Files are names `TR<LOCUS><SEGMENT>.fa`.
#####Receptor type#####
    igblast_seqtype = TCR

Type of sequence to be analysed. Since TraCeR currently only works with TCR sequences, there's no need to change this. 

####Kallisto options####
		[kallisto_options]
		base_transcriptome = /path/to/kallisto/transcriptome

Location of the transcriptome fasta file to which the specific TCR sequences will be appended from each cell. Can be downloaded from http://bio.math.berkeley.edu/kallisto/transcriptomes/ and many other places.

##Using TraCeR##
Tracer has two modes *assemble* and *summarise*. 

*Assemble* takes fastq files of paired-end RNA-seq reads from a single-cell and reconstructs TCR sequences.

*Summarise* takes a set of directories containing output from the *assemble* phase (each directory represents a single cell) and summarises TCR recovery rates as well as generating clonotype networks. 


###*Assemble*: TCR reconstruction##

####Usage####

    ./tracer assemble [options] <file_1> <file_2> <cell_name> <output_directory>

#####Main arguments#####
`<file_1>` : fastq file containing #1 mates from paired-end sequencing  
`<file_2>` : fastq file containing #2 mates from paired-end sequencing   
`<cell_name>` : name of the cell. This is arbitrary text that will be used for all subsequent references to the cell in filenames/labels etc.  
`<output_directory>` : directory for output. Will be created if it doesn't exist. Cell-specific output will go into /<output_directory>/cell_name  

#####Options#####
`-p/--ncores <int>` : number of processor cores available. This is passed to Bowtie2 and Trinity. Default=1.  
`-c/--config_file <conf_file>` : config file to use. Default = `tracer.conf`  
`-r/--resume_with_existing_files` : if this is set, TraCeR will look for existing output files and not re-run steps that already appear to have been completed. This saves time if TraCeR died partway through a step and you want to resume where it left off.   

