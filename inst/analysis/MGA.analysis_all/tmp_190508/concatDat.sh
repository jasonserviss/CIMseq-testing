#!/bin/bash -l
#SBATCH -J concat
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_all/logs/concat.out
#SBATCH -t 00:60:00
#SBATCH -n 1
#SBATCH -A snic2019-3-84
#SBATCH -p devcore
Rscript --vanilla ~/scripts/uppmax.concat.swarm.R /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_all/logs/seedFile_190507.txt
