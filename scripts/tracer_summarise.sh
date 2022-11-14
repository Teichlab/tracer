#!/bin/bash

set -euo pipefail

#Outdir needs to be a singular directory containing sudirectories each containing output for "tracer assemble" for a single sample
#i.e. out-SAMPLE1, out-SAMPLE2, out-SAMPLE3 each containing: "aligned_reads  expression_quantification  filtered_TCR_seqs  IgBLAST_output  Trinity_output  unfiltered_TCR_seqs"

OUTDIR=$1

tracer summarise  --loci A B D G -p 16 -s Hsap -c /home/.tracerrc $OUTDIR
