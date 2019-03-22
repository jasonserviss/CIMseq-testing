#!/bin/bash -l
#SBATCH -J CIMseq-197
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_manualColonSplit/tmp/CIMseq-197-%j.out
#SBATCH -t 72:00:00
#SBATCH -n 1
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_manualColonSplit/scripts/CIMseqSwarm.R 197
