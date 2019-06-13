#!/bin/bash -l
#SBATCH -J concat
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_engeSimpleSynMul/logs/concat.out
#SBATCH -t 00:10:00
#SBATCH -n 1
#SBATCH -A snic2018-8-151
#SBATCH -p devcore
Rscript --vanilla ~/scripts/uppmax.concat.swarm.R /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_engeSimpleSynMul/logs/seedFile_190529.txt
