# TraCeR
TraCeR - reconstruction of T cell receptor sequences from single-cell RNA-seq data.

##Contents##
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Setup](#setup)
4. [Usage](#using-tracer)
	- [*Assemble*](#assemble-tcr-reconstruction)
    - [*Summarise*](#summarise-summary-and-clonotype-networks)


##Introduction
This tool reconstructs the sequences of rearranged and expressed T cell receptor genes from single-cell RNA-seq data. It then uses the TCR sequences to identify cells that have the same receptor sequences and so derive from the same original clonally-expanded cell. 

For more information on TraCeR, its validation and how it can be applied to investigate T cell populations during infection, see our [manuscript on bioRxiv](http://biorxiv.org/content/early/2015/08/28/025676).

Please email questions / problems to mstubb@ebi.ac.uk

##Installation
TraCeR is written in Python and so can just be downloaded, made executable (with `chmod u+x tracer`) and run or run with `python tracer`. Download the latest version and accompanying files from www.github.com/teichlab/tracer. 

Tracer relies on several additional tools and Python modules that you should install.

###Pre-requisites

####Software####
1. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - required for alignment of reads to synthetic TCR genomes.
2. [Trinity](http://sourceforge.net/projects/trinityrnaseq/files/PREV_CONTENTS/previous_releases/) - required for assembly of reads into TCR contigs. **Currently TraCeR uses Trinity parameters intended for use with [Trinity v1](http://sourceforge.net/projects/trinityrnaseq/files/PREV_CONTENTS/previous_releases/). Updates for use with [Trinity v2](https://github.com/trinityrnaseq/trinityrnaseq/wiki) are coming soon.**
3. [IgBLAST](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone) - required for analysis of assembled contigs. [FTP site](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/).
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

**Important:** If you  specify relative paths in the config file these will be used as relative to the directory that contains the `tracer` executable.

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
    trinity_grid_conf = /path/to/trinity/grid.conf

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

Location of the transcriptome fasta file to which the specific TCR sequences will be appended from each cell. Can be downloaded from http://bio.math.berkeley.edu/kallisto/transcriptomes/ and many other places. This must be a plain-text fasta file so decompress it if necessary (files from the Kallisto link are gzipped).

##Using TraCeR##
Tracer has two modes *assemble* and *summarise*. 

*Assemble* takes fastq files of paired-end RNA-seq reads from a single-cell and reconstructs TCR sequences.

*Summarise* takes a set of directories containing output from the *assemble* phase (each directory represents a single cell) and summarises TCR recovery rates as well as generating clonotype networks. 


###*Assemble*: TCR reconstruction##

####Usage####

    tracer assemble [options] <file_1> <file_2> <cell_name> <output_directory>

#####Main arguments#####
`<file_1>` : fastq file containing #1 mates from paired-end sequencing  
`<file_2>` : fastq file containing #2 mates from paired-end sequencing   
`<cell_name>` : name of the cell. This is arbitrary text that will be used for all subsequent references to the cell in filenames/labels etc.  
`<output_directory>` : directory for output. Will be created if it doesn't exist. Cell-specific output will go into `/<output_directory>/<cell_name>`. This path should be the same for every cell that you want to summarise together.

#####Options#####
`-p/--ncores <int>` : number of processor cores available. This is passed to Bowtie2 and Trinity. Default=1.  
`-c/--config_file <conf_file>` : config file to use. Default = `tracer.conf`  
`-s/--species` : Species from which the T cells were derived. Options are `Mmus` or `Hsap` for mouse or human data. This is only important for determination of iNKT cells in the `summarise` step because it defines the V segments that are indicative of iNKT cells. Default = `Mmus`.  
`-r/--resume_with_existing_files` : if this is set, TraCeR will look for existing output files and not re-run steps that already appear to have been completed. This saves time if TraCeR died partway through a step and you want to resume where it left off.   


####Output####

For each cell, an `/<output_directory>/<cell_name>` directory will be created. This will contain the following subdirectories.

1. `<output_directory>/<cell_name>/aligned_reads`  
    This contains the output from Bowtie2 with the sequences of the reads that aligned to the synthetic genomes.

2. `<output_directory>/<cell_name>/Trinity_output`  
    Contains fasta files for each locus where contigs could be assembled. Also two text files that log successful and unsuccessful assemblies.

3. `<output_directory>/<cell_name>/IgBLAST_output`  
    Files with the output from IgBLAST for the contigs from each locus. 

4. `<output_directory>/<cell_name>/unfiltered_TCR_seqs`  
    Files describing the TCR sequences that were assembled prior to filtering by expression if necessary.
    - `unfiltered_TCRs.txt` : text file containing TCR details. Begins with count of productive/total rearrangements detected for each locus. Then details of each detected recombinant.
    - `<cell_name>_TCRseqs.fa` : fasta file containing full-length, reconstructed TCR sequences.
    - `<cell_name>.pkl` : Python [pickle](https://docs.python.org/2/library/pickle.html) file containing the internal representation of the cell and its recombinants as used by TraCeR. This is used in the summarisation steps.

5. `<output_directory>/<cell_name>/expression_quantification`  
    Contains Kallisto output with expression quantification of the entire transcriptome *including* the reconstructed TCRs.

6. `<output_directory>/<cell_name>/filtered_TCR_seqs`  
    Contains the same files as the unfiltered directory above but these recombinants have been filtered so that only the two most highly expressed from each locus are retained. This resolves biologically implausible situtations where more than two recombinants are detected for a locus. **This directory contains the final output with high-confidence TCR assignments**.


###*Summarise*: Summary and clonotype networks###

####Usage####
    tracer summarise [options] <input_dir>

#####Main argument#####
`<input_dir>` : directory containing subdirectories of each cell you want to summarise. 

#####Options#####
`-c/--config_file <conf_file>` : config file to use. Default = `tracer.conf`  
`-u/--use_unfiltered` : Set this flag to use unfiltered recombinants for summary and networks rather than the recombinants filtered by expression level.  
`-i/--keep_inkt` : TraCeR attempts to identify iNKT cells by their characteristic TCRA gene segments (TRAV11–TRAJ18). By default, these are removed before creation of clonotype networks. Setting this option retains the iNKT cells in all stages.    
`-g/--graph_format` : Output format for the clonotype networks. This is passed directly to Graphviz and so must be one of the options detailed at http://www.graphviz.org/doc/info/output.html. 

####Output####
Output is written to `<input_dir>/filtered_TCR_summary` or `<input_dir>/unfiltered_TCR_summary` depending on whether the `--use_unfiltered` option was set.

The following output files are generated:

1. `TCR_summary.txt`
    Summary statistics describing successful TCR reconstruction rates and the numbers of cells with 0, 1, 2 or more recombinants for each locus.
2. `reconstructed_lengths_TCR[A|B].pdf`
    Distribution plots showing the lengths of the VDJ regions from assembled TCR contigs. Longer contigs give higher-confidence segment assignments.
3. `clonotype_sizes.pdf`
    Distribution of clonotype sizes.
4.  `clonotype_network_[with|without]_identifiers.<graph_format>`
    graphical representation of clonotype networks either with full recombinant identifiers or just lines indicating presence/absence of recombinants.
5.  `clonotype_network_[with|without]_identifiers.dot`
    files describing the clonotype networks in the [Graphviz DOT language](http://www.graphviz.org/doc/info/lang.html)
