#!/bin/bash -l
#SBATCH -J concat
#SBATCH -o /proj/snic2018-8-151/private/Github/CIMseq.testing/inst/analysis/MGA.analysis_SI_strictClasses/logs/concat.out
#SBATCH -t 00:60:00
#SBATCH -n 1
#SBATCH -A snic2018-8-151
#SBATCH -p devcore
Rscript --vanilla ~/scripts/uppmax.concat.swarm.R 435
