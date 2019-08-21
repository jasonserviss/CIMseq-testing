#!/bin/bash -l
#SBATCH -J CIMseq-s.NJA01803.O20
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/deconvoluteSingletsMouse/seedFile_190730/CIMseq-s.NJA01803.O20-%j.out
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/deconvoluteSingletsMouse/scripts/CIMseqSwarm.R s.NJA01803.O20 seedFile_190730
