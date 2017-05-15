# TraCeR
TraCeR - reconstruction of T cell receptor sequences from single-cell RNA-seq data.

**IMPORTANT: Python dependencies have changed since the last release. Use the requirements file (detailed [here](#setup)) to update them if necessary.**

## Contents ##
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Setup](#setup)
4. [Testing](#testing-tracer)
5. [Usage](#using-tracer)
	- [*Assemble*](#assemble-tcr-reconstruction)
    - [*Summarise*](#summarise-summary-and-clonotype-networks)


## Introduction
This tool reconstructs the sequences of rearranged and expressed T cell receptor genes from single-cell RNA-seq data. It then uses the TCR sequences to identify cells that have the same receptor sequences and so derive from the same original clonally-expanded cell. 

For more information on TraCeR, its validation and how it can be applied to investigate T cell populations during infection, see our [paper in Nature Methods](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3800.html) or the [bioRxiv preprint](http://biorxiv.org/content/early/2015/08/28/025676) that preceded it.

Please email questions / problems to ms31@sanger.ac.uk

## Installation
TraCeR is written in Python and so can just be downloaded, made executable (with `chmod u+x tracer`) and run or run with `python tracer`. Download the latest version and accompanying files from www.github.com/teichlab/tracer. 

TraCeR relies on several additional tools and Python modules that you should install.

Note that TraCeR is compatible with both Python 2 and 3.

### Pre-requisites

#### Software 
1. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - required for alignment of reads to synthetic TCR genomes.
2. [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) - required for assembly of reads into TCR contigs. TraCeR now works with both version 1 and version 2 of Trinity. It should automatically detect the version that is installed or you can [specify it in the config file](https://github.com/Teichlab/tracer#trinity-options).
    - Please note that Trinity requires a working installation of [Bowtie v1](http://bowtie-bio.sourceforge.net).
3. [IgBLAST](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone) - required for analysis of assembled contigs. (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/).
4. [makeblastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ ) - **optional** but required if you want to use TraCeR's `build` mode to make your own references.
5. Software for quantification of TCR expression:
    * [Kallisto](http://pachterlab.github.io/kallisto/), or alternatively
    * [Salmon](https://github.com/COMBINE-lab/salmon/releases).
6. [Graphviz](http://www.graphviz.org) - Dot and Neato drawing programs required for visualisation of clonotype graphs. This is optional - see the [`--no_networks` option](#options-1) to [`summarise`](#summarise-summary-and-clonotype-networks).

##### Installing IgBlast 
Downloading the executable files from `ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/<version_number>` is not sufficient for a working IgBlast installation. You must also download the `internal_data` directory (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data) and put it into the same directory as the igblast executable. This is also described in the igblast README file.

You should also ensure to set the `$IGDATA` environment variable to point to the location of the IgBlast executable. For example run `export IGDATA=/<path_to_igblast>/igblast/1.4.0/bin`.

#### Python modules 
1. [Matplotlib](http://matplotlib.org)
2. [Seaborn](http://stanford.edu/~mwaskom/software/seaborn/)
3. [Biopython](http://biopython.org/)
4. [Prettytable](https://code.google.com/p/prettytable/)
5. [Levenshtein](https://pypi.python.org/pypi/python-Levenshtein/)
6. [Networkx](https://networkx.github.io)
7. It seems that v1.11 of Networkx behaves differently when writing dot files for use with Graphviz. If you have this (or later versions) you also need to install [PyDotPlus](http://pydotplus.readthedocs.org).
8. [Future](http://python-future.org/index.html) for compatibility with Python 2.

## Setup 
To set up the python dependencies, use the requirements file:

    pip install -r requirements.txt

It is **highly** recommended that numpy and biopython are first installed through your system's package manager or conda.

The tracer module is then installed using:

    python setup.py install

This will add the binary 'tracer' to your local bin folder, which can then be run from anywhere.

If you would like to contribute to TraCeR, you can set up a development version with

    python setup.py develop

Which will make TraCeR accessible in your python environment, and incorporate local updates to the code.

Once the prerequisites above are installed and working you're ready to tell TraCeR where to find them.

TraCeR uses a configuration file to point it to the locations of files that it needs and a couple of other options.
An example configuration file is included in the repository - `tracer.conf`.
By default, this is `~/.tracerrc`. If tracer fails to find this file, it will use the `tracer.conf` in the repository.
 The `-c` option to the various tracer modules allows you to specify any other file to act as the configuration file.

**Important:** If you  specify relative paths in the config file these will be used as relative to the main installation directory. For example, `resources/Mmus/igblast_dbs` will resolve to `/<wherever you installed tracer>/tracer/resources/Mmus/igblast_dbs`.

### External tool locations 
Tracer will look in your system's `PATH` for external tools. You can override this behaviour by editing your `~/.tracerrc`.
Edit `~/.tracerrc` (or a copy) so that the paths within the `[tool_locations]` section point to the executables for all of the required tools.

	[tool_locations]
	#paths to tools used by TraCeR for alignment, quantitation, etc
	bowtie2_path = /path/to/bowtie2
	bowtie2-build_path = /path/to/bowtie2-build
	igblast_path = /path/to/igblastn
	makeblastdb_path = /path/to/makeblastdb
	kallisto_path = /path/to/kallisto
	salmon_path = /path/to/salmon
	trinity_path = /path/to/trinity
	dot_path = /path/to/dot
	neato_path = /path/to/neato
		
		
### Resource locations and necessary files 
The tools used by TraCeR need a variety of additional files to work properly and to allow extraction of TCR-derived reads and expression quantification etc. The locations of these files are specified in the other sections of the configuration file and are detailed below.

Currently, organism-specific files (TCR gene sequences, synthetic genome indices, igblast_indices) for mouse and human are distributed with the source-code in the `resources` directory. There will soon be a `build` module that constructs all of the necessary synthetic genomes and indices from any collection of V and J gene sequences. 


#### Bowtie synthetic genomes path 
	[bowtie2_options]
	synthetic_genome_index_path = resources/Mmus/synthetic_genomes

This path specifies the directory that contains Bowtie2 indices constructed from all possible combinations of V and J segments for each locus. 

#### Trinity options 
##### Jellyfish memory 
	[trinity_options]
	#line below specifies maximum memory for Trinity Jellyfish component. Set it appropriately for your environment.
	max_jellyfish_memory = 1G

Trinity needs to know the maximum memory available to it for the Jellyfish component. Specify this here.


#### Trinity version 
    #uncomment the line below to explicitly specify Trinity version. Options are '1' or '2'
    #trinity_version = 2

TraCeR will automatically detect the version of Trinity you have installed. You can also explicitly specify it here if you wish.

##### HPC configuration 
    #uncomment the line below if you've got a configuration file for Trinity to use a computing grid #
    trinity_grid_conf = /path/to/trinity/grid.conf

Trinity can parallelise contig assembly by submitting jobs across a compute cluster. If you're running in such an environment you can specify an optional trinity config file here. See the Trinity documentation for more information.
 

#### IgBLAST options 
##### Databases path 
	[IgBlast_options]
	igblast_index_location = resources/Mmus/igblast_dbs
##### VDJ sequences 
This path specifies the directory that contains IgBLAST database files for V, D and J genes. These files are named `imgt_tcr_db_<SEGMENT>.fa`.

    imgt_seq_location = resources/Mmus/imgt_sequences
		
Path to fasta files with sequences for each V, D or J gene. Files are names `TR<LOCUS><SEGMENT>.fa`.
##### Receptor type 
    igblast_seqtype = TCR

Type of sequence to be analysed. Since TraCeR currently only works with TCR sequences, there's no need to change this. 

#### Base transcriptomes for Kallisto/Salmon 
	[base_transcriptomes]
	Mmus = /path/to/kallisto/transcriptome_for_Mmus
	Hsap = /path/to/kallisto/transcriptome_for_Hsap

Location of the transcriptome fasta file to which the specific TCR sequences will be appended from each cell. Can be downloaded from http://bio.math.berkeley.edu/kallisto/transcriptomes/ and many other places. This must be a plain-text fasta file so decompress it if necessary (files from the Kallisto link are gzipped).

#### Base indices for Kallisto/Salmon
	[salmon_base_indices]
	Mmus = /path/to/salmon/index_for_Mmus
	Hsap = /path/to/salmon/index_for_Hsap
	[kallisto_base_indices]
	Mmus = /path/to/kallisto/index_for_Mmus
	Hsap = /path/to/kallisto/index_for_Hsap

Location of Kallisto/Salmon indices built (exclusively) from corresponding `[base_transcriptomes]`. These indices are only needed when option `--small_index` is used in *Assemble* mode (see below). 

#### Salmon options
	[salmon_options]
	libType = A
	kmerLen = 31

* Description of the type of sequencing library from which the reads come (containing, e.g., the relative orientation of paired end reads). As of version 0.7.0, Salmon also has the ability to automatically infer (i.e. guess) the library type based on how the first few thousand reads map to the transcriptome. Set `libType = A` for automatic detection.
* Salmon builds the quasi-mapping-based index, using an auxiliary k-mer hash over k-mers of length `kmerLen`. While quasi-mapping will make used of arbitrarily long matches between the query and reference, the k size selected here will act as the minimum acceptable length for a valid match. The value for `kmerLen` must be odd; its default and maximum value is 31. 

See salmon [documentation](http://salmon.readthedocs.io/en/latest/salmon.html) for more details.


## Testing TraCeR 
TraCeR comes with a small dataset in `test_data/` (containing only TCRA or TCRB reads for a single cell) that you can use to test your installation and config file and confirm that all the prerequisites are working. Run it as:

    tracer test -p <ncores> -c <config_file>
    
**Note:** The data used in the test are derived from mouse T cells so make sure that the config file points to the appropriate mouse resource files.

You can also pass the following two options to change the Graphviz output format or to prevent attempts to draw network graphs

`-g/--graph_format` : Output format for the clonotype networks. This is passed directly to Graphviz and so must be one of the options detailed at http://www.graphviz.org/doc/info/output.html.  
`--no_networks` : Don't try to draw clonotype network graphs. This is useful if you don't have a working installation of Graphviz.
    
Running `test` will peform the [`assemble`](#assemble-tcr-reconstruction) step using the small test dataset. It will then perform [`summarise`](#summarise-summary-and-clonotype-networks) using the assemblies that are generated along with pre-calculated output for two other cells (in `test_data/results`).

Compare the output in `test_data/results/filtered_TCR_summary` with the expected results in `test_data/expected_summary`. There should be three cells, each with one productive alpha, one productive beta, one non-productive alpha and one non-productive beta. Cells 1 and 2 should be in a clonotype.


## Using TraCeR 
Tracer has three modes: *assemble*, *summarise* and *build*. 

*Assemble* takes fastq files of paired-end RNA-seq reads from a single-cell and reconstructs TCR sequences.

*Summarise* takes a set of directories containing output from the *assemble* phase (each directory represents a single cell) and summarises TCR recovery rates as well as generating clonotype networks. 

*Build* creates new combinatorial recombinomes for species other than the inbuilt Human and Mouse.


### *Assemble*: TCR reconstruction 

#### Usage 

    tracer assemble [options] <file_1> [<file_2>] <cell_name> <output_directory>


##### Main arguments
* `<file_1>` : fastq file containing #1 mates from paired-end sequencing or all reads from single-end sequencing.   
* `<file_2>` : fastq file containing #2 mates from paired-end sequencing. Do not use if your data are from single-end sequencing.  
* `<cell_name>` : name of the cell. This is arbitrary text that will be used for all subsequent references to the cell in filenames/labels etc.     
* `<output_directory>` : directory for output. Will be created if it doesn't exist. Cell-specific output will go into `/<output_directory>/<cell_name>`. This path should be the same for every cell that you want to summarise together.

##### Options 

* `-p/--ncores <int>` : number of processor cores available. This is passed to Bowtie2, Trinity, and Kallisto or Salmon. Default=1.
* `-c/--config_file <conf_file>` : config file to use. Default = `~/.tracerrc`
* `-s/--species` : Species from which the T cells were derived. Options are `Mmus` or `Hsap` for mouse or human data. This is only important for determination of iNKT cells in the `summarise` step because it defines the V segments that are indicative of iNKT cells. Default = `Mmus`.
* `-r/--resume_with_existing_files` : if this is set, TraCeR will look for existing output files and not re-run steps that already appear to have been completed. This saves time if TraCeR died partway through a step and you want to resume where it left off.
* `-m/--seq_method` : method by which to generate sequences for assessment of recombinant productivity. By default (`-m imgt`), TraCeR replaces all but the junctional sequence of each detected recombinant with the reference sequence from IMGT prior to assessing productivity of the sequence. This makes the assumption that sequence changes outside the junctional region are due to PCR/sequencing errors rather than being genuine polymorphisms. This is likely to be true for well-characterised mouse sequences but may be less so for human and other outbred populations. To determine productivity from only the assembled contig sequence for each recombinant use `-m assembly`.
* `-q/--quant_method` : Method used for expression quantification. Options are `-q salmon` and `-q kallisto` (default).
* `--small_index` : Use this option to speed up expression quantification. The location of an index that is built (exclusively) from the corresponding `base_transcriptome` must be specified in the configuration file (under `[salmon_base_indices]` or `[kallisto_base_indices]`). Since this index does not contain any TCR sequences, it has to be built only once for each species and can then be used for all cells. When option `--small-index` is used, reads are first quantified with this `base_index` (this should lead to non-zero TPM numbers for all non-TCR transcripts present in the cell). After selecting all transcripts with non-zero TPM numbers, cell-specific TCR sequences (as constructed by Bowtie-Trinity-IgBlast) are appended to that list. Then a new Kallisto/Salmon index is built from the combined set (since this set contains only a subset of the `base_transcriptome`, this is fast). Finally, reads are quantified with this small index. 
* `--single_end` : use this option if your data are single-end reads. If this option is set you must specify fragment length and fragment sd as below.
* `--fragment_length` : Estimated average fragment length in the sequencing library. Used for Kallisto quantification. Required for single-end data. Can also be set for paired-end data if you don't want Kallisto to estimate it directly.
* `--fragment_sd` : Estimated standard deviation of average fragment length in the sequencing library. Used for Kallisto quantification. Required for single-end data. Can also be set for paired-end data if you don't want Kallisto to estimate it directly.
* `--invariant_sequences`: Custom invariant sequence file. Use the default example in 'resources/Mmus/invariant_seqs.csv'

#### Output 

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
    Contains Kallisto/Salmon output with expression quantification of the entire transcriptome *including* the reconstructed TCRs. When option `--small_index` is used, this directory contains only the output of the quantification with the small index (built from reconstructed TCRs and only a subset of the base transcriptome; see above).

6. `<output_directory>/<cell_name>/filtered_TCR_seqs`  
    Contains the same files as the unfiltered directory above but these recombinants have been filtered so that only the two most highly expressed from each locus are retained. This resolves biologically implausible situtations where more than two recombinants are detected for a locus. **This directory contains the final output with high-confidence TCR assignments**.


### *Summarise*: Summary and clonotype networks 

#### Usage 
    tracer summarise [options] <input_dir>

##### Main argument 
* `<input_dir>` : directory containing subdirectories of each cell you want to summarise. 

##### Options 
* `-c/--config_file <conf_file>` : config file to use. Default = `~/.tracerrc`
* `-u/--use_unfiltered` : Set this flag to use unfiltered recombinants for summary and networks rather than the recombinants filtered by expression level.  
* `-i/--keep_inkt` : TraCeR attempts to identify iNKT cells by their characteristic TCRA gene segments (TRAV11–TRAJ18). By default, these are removed before creation of clonotype networks. Setting this option retains the iNKT cells in all stages.    
* `-g/--graph_format` : Output format for the clonotype networks. This is passed directly to Graphviz and so must be one of the options detailed at http://www.graphviz.org/doc/info/output.html.  
* `--no_networks` : Don't try to draw clonotype network graphs. This is useful if you don't have a working installation of Graphviz.

#### Output 
Output is written to `<input_dir>/filtered_TCR_summary` or `<input_dir>/unfiltered_TCR_summary` depending on whether the `--use_unfiltered` option was set.

The following output files are generated:

1. `TCR_summary.txt`
    Summary statistics describing successful TCR reconstruction rates and the numbers of cells with 0, 1, 2 or more recombinants for each locus.
2. `recombinants.txt`
    List of TCR identifiers, lengths and productivities for each cell. 
3. `reconstructed_lengths_TCR[A|B].pdf` and  `reconstructed_lengths_TCR[A|B].txt`
    Distribution plots (and text files with underlying data) showing the lengths of the VDJ regions from assembled TCR contigs. Longer contigs give higher-confidence segment assignments. Text files are only generated if at least one TCR is found for a locus. Plots are only generated if at least two TCRs are found for a locus. 
4. `clonotype_sizes.pdf` and `clonotype_sizes.txt`
    Distribution of clonotype sizes as bar graph and text file.
5.  `clonotype_network_[with|without]_identifiers.<graph_format>`
    graphical representation of clonotype networks either with full recombinant identifiers or just lines indicating presence/absence of recombinants.
6.  `clonotype_network_[with|without]_identifiers.dot`
    files describing the clonotype networks in the [Graphviz DOT language](http://www.graphviz.org/doc/info/lang.html)
 
### *Build*: Build Combinatorial Recombinomes for a Given Species

#### Usage
    tracer build <species> <receptor_name> <locus_name> <N_padding> <colour> <V_seqs> <J_seqs> <C_seqs> <D_seqs>

#### Main Arguments
* `<species>` : Species (e.g. Mmus).
* `<receptor_name>` : Name of receptor (e.g. TCR).
* `<locus_name>` : Name of locus (e.g. A)
* `<N_padding>` : Number of ambiguous N nucleotides between V and J
* `<colour>` : Colour for the productive recombinants (optional). Specify as HTML (e.g. E41A1C) or use "random"
* `<V_seqs>` : Fasta file containing V gene sequences
* `<J_seqs>` : Fasta file containing J gene sequences
* `<C_seqs>` : Fasta file containing single constant region sequence
* `<D_seqs>` : Fasta file containing D gene sequences (optional)

#### Options

* `-f/--force_overwrite` : Force overwrite of existing resources
