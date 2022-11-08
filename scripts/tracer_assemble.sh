#!/bin/bash

set -euo pipefail

##SANGER EXAMPLE: 
##image=/nfs/cellgeni/singularity/images/tracer-ubuntu-20_04-trinityrnaseq-2_14_0-igblast-1_19_0-kallisto-0_48_0-grch38-transcript_41-grcm38-transcript_m30.sif

image=/path/to/tracer/image.sif

#fastq directory containing subdirectories of each sample, each sample subdirectory contains a single R1 fastq and a single R2 fastq
FQDIR=$1

#sample file containing list of samples that correspond to sample directories within FQDIR
SAMPLE_FILE=$2

cat $SAMPLE_FILE | while read SAMPLE; do
  mkdir -p "tracer_output_full/out-${SAMPLE}"
  for fq in "${FQDIR}/${SAMPLE}/"*"fastq.gz"; do 
    if [[ "$fq" == *"R1"* ]]; then
      FQ1=$fq
    elif [[ "$fq" == *"R2"* ]]; then
      FQ2=$fq
    fi
  done
  zcat ${FQ1} > ${SAMPLE}.R1.fastq 
  zcat ${FQ2} > ${SAMPLE}.R2.fastq
  #Please note if you want to mount multiple paths from the system youre running tracer on to the singularity image then you need to do -B /path1/to/mount,/path2/to/mount
  singularity exec -B /path/to/mount ${image} tracer assemble --loci A B D G -p 4 -s Hsap -c /home/.tracerrc ${FQ1} ${FQ2} out-${SAMPLE} /path/to/output/directory
  ##SANGER EXAMPLE: 
  ##singularity exec -B /lustre,/nfs ${image} tracer assemble --loci A B D G -p 4 -s Hsap -c /home/.tracerrc ${FQ1} ${FQ2} out-${SAMPLE} tracer_output
  rm ${SAMPLE}.R1.fastq
  rm ${SAMPLE}.R2.fastq
done
