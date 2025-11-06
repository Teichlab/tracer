# TraCeR
TraCeR - reconstruction of T cell receptor sequences from single-cell RNA-seq data.

## Contents ##
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Setup](#setup)
4. [Testing](#testing-tracer)
5. [Usage](#using-tracer)
	- [*Assemble*](#assemble-tcr-reconstruction)
    - [*Summarise*](#summarise-summary-and-clonotype-networks)
6. [Docker image](#docker-image)
7. [Wrapper](#wrapper)


## Introduction
This tool reconstructs the sequences of rearranged and expressed T cell receptor genes from single-cell RNA-seq data. It then uses the TCR sequences to identify cells that have the same receptor sequences and so derive from the same original clonally-expanded cell. 

For more information on TraCeR, its validation and how it can be applied to investigate T cell populations during infection, see our [paper in Nature Methods](http://dx.doi.org/10.1038/nmeth.3800) or the [bioRxiv preprint](http://biorxiv.org/content/early/2015/08/28/025676) that preceded it.

Please email questions / problems to mjt.stubbington \[at\] gmail.com.

_Also, please note that I (Mike) created and developed TraCeR while I was at EMBL-EBI and subsequently at the Sanger Institute. Now that I have left Sanger, I will continue to maintain the code when possible **as a personal project** independent of any employment that I have_

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
Downloading the executable files from `ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/<version_number>` is not sufficient for a working IgBlast installation. You must also download the `internal_data` directory (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data) and `optional_file` directory (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/optional_file/) and put them into the same directory as the igblast executable. This is also described in the igblast README file.

You should also ensure to set the `$IGDATA` environment variable to point to the location of the IgBlast executable. For example run `export IGDATA=/<path_to_igblast>/igblast/1.4.0/bin`.

## Setup 
To set up the python dependencies, use the requirements file:

    pip install -r requirements.txt

It is **highly** recommended that numpy and biopython are first installed through your system's package manager or conda.

Note: Seaborn depends on the module statsmodels, which if updated through other packages may cause problems in Seaborn. If such issues arise, try to uninstall statsmodels and install again:
  
    conda uninstall statsmodels --yes
    conda install -c taugspurger statsmodels=0.8.0     

The tracer module is then installed using:

    python setup.py install

This will add the binary 'tracer' to your local bin folder, which can then be run from anywhere.

_Note that installing TraCeR using this method requires you to specify the location of the originally downloaded files in your config file ([see below](#tracer-directory))._

If you would like to contribute to TraCeR, you can set up a development version with

    python setup.py develop

Which will make TraCeR accessible in your python environment, and incorporate local updates to the code.

Once the prerequisites above are installed and working you're ready to tell TraCeR where to find them.

TraCeR uses a configuration file to point it to the locations of files that it needs and a couple of other options.
An example configuration file is included in the repository - `tracer.conf`.

TraCeR looks for the configuration file, in descending order of priority, from the following sources:
1. The `-c` option used at run time for any of the TraCeR modes
2. The environmental variable `TRACER_CONF`, which can be set to point to a default location
3. The default global location of `~/.tracerrc`
4. If all else fails, the provided example `tracer.conf` in the directory is used

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

#### Trinity with short reads (<50 bases)

If you're running Trinity with read lengths that are shorter than 50 bases, you'll be restricted to using the Inchworm component of Trinity, which does draft contig assembly via greedy kmer extension. In your 'tracer.conf', uncomment the 'inchworm_only=True' to activate short read mode, and uncomment the 'trinity_kmer_length' as well.  You can change the trinity_kmer_length value as needed; 17 is the shortest value, and the default length in Trinity is 25.  Trinity v2 should be used in this mode.



##### HPC configuration 
    #uncomment the line below if you've got a configuration file for Trinity to use a computing grid #
    trinity_grid_conf = /path/to/trinity/grid.conf

Trinity can parallelise contig assembly by submitting jobs across a compute cluster. If you're running in such an environment you can specify an optional trinity config file here. See the Trinity documentation for more information.
 

#### IgBLAST options 
##### Receptor type 
    igblast_seqtype = TCR

Type of sequence to be analysed. Since TraCeR currently only works with TCR sequences, there's no need to change this. 

#### Base transcriptomes for Kallisto/Salmon 
	[base_transcriptomes]
	Mmus = /path/to/kallisto/transcriptome_for_Mmus
	Hsap = /path/to/kallisto/transcriptome_for_Hsap

Location of the transcriptome fasta file to which the specific TCR sequences will be appended from each cell. Instructions on how to make reference transcriptomes can be found at: http://www.nxn.se/valent/2016/10/3/the-first-steps-in-rna-seq-expression-analysis-single-cell-and-other. Reference trasncriptomes should be uncompressed fasta files.
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

### TraCeR directory

        [tracer_location]
        #Path to where TraCeR was originally downloaded
        tracer_path = /path/to/tracer

Location of the cloned TraCeR repository containing tracerlib, test_data, resources etc. Eg. `/user/software/tracer`.
Needed for localisation of resources and test_data if running TraCeR with the tracer binary. Note: this is the path to the **directory** and not the executable.

## Testing TraCeR 
TraCeR comes with a small dataset in `test_data/` (containing only TCRA or TCRB reads for a single cell) that you can use to test your installation and config file and confirm that all the prerequisites are working. Run it as:

    tracer test -p <ncores> -c <config_file>
    
**Note:** The data used in the test are derived from mouse T cells so make sure that the config file points to the appropriate mouse resource files.

You can also pass the following two options to change the Graphviz output format or to prevent attempts to draw network graphs

`-g/--graph_format` : Output format for the clonotype networks. This is passed directly to Graphviz and so must be one of the options detailed at http://www.graphviz.org/doc/info/output.html.  
`--no_networks` : Don't try to draw clonotype network graphs. This is useful if you don't have a working installation of Graphviz.
    
Running `test` will peform the [`assemble`](#assemble-tcr-reconstruction) step using the small test dataset.
 It will then perform [`summarise`](#summarise-summary-and-clonotype-networks) using the assemblies that are generated along with pre-calculated output for two other cells (in `test_data/results`), dumping the data in the `output` directory as specified.

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
* `--resource_dir`: the directory containing the resources required for alignment. By default this is the `resources` directory in this repository, but can be pointed to a user-built set of resources.
* `-s/--species` : Species from which the T cells were derived. Options are `Mmus` or `Hsap` for mouse or human data. If you have defined new species using the [`build`](#build-build-combinatorial-recombinomes-for-a-given-species) mode, you should specify the same name here. Default = `Mmus`.
* `--receptor_name` : If, for some reason, you've used a different receptor name when using [`build`](#build-build-combinatorial-recombinomes-for-a-given-species) then you should also specify it here. Default = `TCR`
* `--loci` : The specific loci that you wish to assemble. TraCeR knows about alpha (`A`), beta (`B`), gamma (`G`) and delta (`D`) for mouse and human TCRs. Other locus names can be specified if you use [`build`](#build-build-combinatorial-recombinomes-for-a-given-species). By default, TraCeR will attempt to assemble alpha and beta sequences. To change this, pass a space-delimited list of locus names. For example, to only look for gamma/delta sequences use `--loci G D`. **Important** The way that this argument is parsed means that you **must** follow it with another command line option or the argument parser will get confused. For example: `--loci G D --p 1`. We know this [isn't ideal](https://github.com/Teichlab/tracer/issues/62).
* `-r/--resume_with_existing_files` : if this is set, TraCeR will look for existing output files and not re-run steps that already appear to have been completed. This saves time if TraCeR died partway through a step and you want to resume where it left off.
* `-m/--seq_method` : method by which to generate sequences for assessment of recombinant productivity. By default (`-m imgt`), TraCeR replaces all but the junctional sequence of each detected recombinant with the reference sequence from IMGT prior to assessing productivity of the sequence. This makes the assumption that sequence changes outside the junctional region are due to PCR/sequencing errors rather than being genuine polymorphisms. This is likely to be true for well-characterised mouse sequences but may be less so for human and other outbred populations. To determine productivity from only the assembled contig sequence for each recombinant use `-m assembly`.
* `-q/--quant_method` : Method used for expression quantification. Options are `-q salmon` and `-q kallisto` (default).
* `--small_index` : Use this option to speed up expression quantification. The location of an index that is built (exclusively) from the corresponding `base_transcriptome` must be specified in the configuration file (under `[salmon_base_indices]` or `[kallisto_base_indices]`). Since this index does not contain any TCR sequences, it has to be built only once for each species and can then be used for all cells. When option `--small-index` is used, reads are first quantified with this `base_index` (this should lead to non-zero TPM numbers for all non-TCR transcripts present in the cell). After selecting all transcripts with non-zero TPM numbers, cell-specific TCR sequences (as constructed by Bowtie-Trinity-IgBlast) are appended to that list. Then a new Kallisto/Salmon index is built from the combined set (since this set contains only a subset of the `base_transcriptome`, this is fast). Finally, reads are quantified with this small index. 
* `--single_end` : use this option if your data are single-end reads. If this option is set you must specify fragment length and fragment sd as below.
* `--fragment_length` : Estimated average fragment length in the sequencing library. Used for Kallisto quantification. Required for single-end data. Can also be set for paired-end data if you don't want Kallisto to estimate it directly.
* `--fragment_sd` : Estimated standard deviation of average fragment length in the sequencing library. Used for Kallisto quantification. Required for single-end data. Can also be set for paired-end data if you don't want Kallisto to estimate it directly.
* `--max_junc_len` : Maximum permitted length of junction string in recombinant identifier. Used to filter out artefacts. May need to be longer for TCR Delta.
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


#### Usage 
    tracer summarise [options] <input_dir>

##### Main argument 
* `<input_dir>` : directory containing subdirectories of each cell you want to summarise. 
* `<resource_dir>`: the directory containing the resources required for alignment.
 By default this is the `resources` directory in this repository, but can be pointed to a user-built set of resources.

##### Options 
* `-c/--config_file <conf_file>` : config file to use. Default = `~/.tracerrc`
* `-u/--use_unfiltered` : Set this flag to use unfiltered recombinants for summary and networks rather than the recombinants filtered by expression level.  
* `--resource_dir`: the directory containing the resources required for alignment. By default this is the `resources` directory in this repository, but can be pointed to a user-built set of resources.
* `-s/--species` : Species from which the T cells were derived. Options are `Mmus` or `Hsap` for mouse or human data. If you have defined new species using the [`build`](#build-build-combinatorial-recombinomes-for-a-given-species) mode, you should specify the same name here. Default = `Mmus`.
* `--receptor_name` : If, for some reason, you've used a different receptor name when using [`build`](#build-build-combinatorial-recombinomes-for-a-given-species) then you should also specify it here. Default = `TCR`
* `--loci` : The specific loci that you wish to summarise. TraCeR knows about alpha (`A`), beta (`B`), gamma (`G`) and delta (`D`) for mouse and human TCRs. Other locus names can be specified if you use [`build`](#build-build-combinatorial-recombinomes-for-a-given-species). By default, TraCeR will attempt to summarise alpha and beta sequences. To change this, pass a space-delimited list of locus names. For example, to only look for gamma/delta sequences use `--loci G D`. **Important** The way that this argument is parsed means that you **must** follow it with another command line option or the argument parser will get confused. For example: `--loci G D --p 1`. We know this [isn't ideal](https://github.com/Teichlab/tracer/issues/62). Default = `A B`
* `-i/--keep_invariant` : TraCeR attempts to identify invariant TCR cells by their characteristic TCRA gene segments. By default, these are removed before creation of clonotype networks. Setting this option retains the invariant cells in all stages.    
* `-g/--graph_format` : Output format for the clonotype networks. This is passed directly to Graphviz and so must be one of the options detailed at http://www.graphviz.org/doc/info/output.html.  
* `--no_networks` : Don't try to draw clonotype network graphs. This is useful if you don't have a working installation of Graphviz.

#### Output 
Output is written to `<input_dir>/filtered_TCR_summary` or `<input_dir>/unfiltered_TCR_summary` depending on whether the `--use_unfiltered` option was set.

The following output files are generated:

1. `TCR_summary.txt`
    Summary statistics describing successful TCR reconstruction rates and the numbers of cells with 0, 1, 2 or more recombinants for each locus.
2. `recombinants.txt`
    List of TCR identifiers, lengths, productivities, and CDR3 sequences (nucleotide and amino acid) for each cell. **Note:** It's possible for non-productive rearrangements to still have a detectable CDR3 if a frameshift introduces a stop-codon after this region. In these cases, the CDR3 nucleotides and amino acids are still reported but are shown in lower-case.
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
* `-o/--output_dir` : Optional directory to write new resource files to. If not specified, the built resources will be placed in the default TraCeR resource directory.  

## Docker image

TraCeR is also available as a standalone Docker image on [DockerHub](https://hub.docker.com/r/teichlab/tracer/), with all of its dependencies installed and configured appropriately. Running TraCeR from the image is very similar to running it from a normal installation. You can pass all the appropriate arguments to the Docker commandwith the usual syntax as described above. Two differences from installing it yourself are that:

* `--small_index` is not supported
* You don't need to worry about specifying a configuration file. This is included in the container.

To run the TraCeR Docker image, run the following command from within a directory that contains your input data:

	docker run -it --rm -v $PWD:/scratch -w /scratch teichlab/tracer
	
followed by any arguments that you want to pass to TraCeR. 

The `-it` flag ensures that you see all the information TraCeR prints to the screen during its run, while `--rm` deletes the individual container created based on the image for the purpose of the run once the analysis is complete, not littering your computer's drive. `-v` creates a volume, allowing the created container to see the contents of your current working directory, and the `-w` flag sets the container's working directory to the newly created volume. 

For example, if you wanted to run the test analysis, you should clone this GitHub repository, navigate to its main directory so you can see the `test_data` folder, and call the following (you need to specify the `-o test_data` so that the results get written to the volume you created, ensuring you can see them after the analysis is finished):

	docker run -it --rm -v $PWD:/scratch -w /scratch teichlab/tracer test -o test_data

If you wish to use `tracer build`, you will need to specify `--output_dir /scratch`, as otherwise the resulting resources will be saved in the default location of the container and subsequently get forgotten about when the build analysis completes, making them unuseable for any actual analyses you may want to perform. This will make the Docker container save the resulting resources in the volume you created, and you can use them for assemble/summarise by running the Dockerised TraCeR from the same directory as the one you used for the build and specifying `--resource_dir /scratch`.

You may need to explicitly tell Docker to increase the memory that it can use. Instructions for [Windows](https://docs.docker.com/docker-for-windows/#advanced) and [Mac](https://docs.docker.com/docker-for-mac/#advanced). Something like 6 or 8 GB is likely to be ok.

## Wrapper

To accelerate the analysis of SmartSeq2 datasets, where each well in a 384 well-plate corresponds to a cell, we have developed a wrapper function, parallelising the `tracer assemble` function.

#### Usage 
    ./scripts/tracer_SmartSeq2.sh [options] <fastq file directory>

##### Main argument 
* `<fastq file directory>` : directory containing all the fastq files from a SmartSeq2 experiment

##### Options 
* `-f` : config file to use. Default: Docker container specific configuration file  
* `-r`: the directory containing the resources required for alignment. By default this is the `resources` directory in this repository, but can be pointed to a user-built set of resources.
* `-a` : Species from which the T cells were derived. Options are `Mmus` or `Hsap` for mouse or human data. If you have defined new species using the [`build`](#build-build-combinatorial-recombinomes-for-a-given-species) mode, you should specify the same name here. Default = `Mmus`.
* `-c` :number of cores    : Number of cores, which TRACER should use. Default = 1.
* `-i` :identifier in the fastq name, distinguishing Read 1 and Read 2. Default = "_R".
* `-s` :single end?: use this flag, if single end sequencing strategy is used. Default = false.
* `-l` :if single end is used: determine fragment length. Default = 75.
* `-d` :if single end is used: determine fragment standard deviation: Default 2.
#

#### Output 
If multiple lanes were used to generate the fastq they will be merged und unzipped. These processed fastq files will be placed to `<outdir>/final`. Final Output is written to `<outdir>/final`. a folder is generated for each cell and the contents are the same as described for the assemble function.

