#!/bin/bash -l
#SBATCH -J CIMseq-36
#SBATCH -o /proj/snic2018-8-151/private/Github/CIMseq.testing/inst/analysis/MGA.analysis_SI_strictClasses/tmp/CIMseq-36-%j.out
#SBATCH -t 18:00:00
#SBATCH -n 1
#SBATCH -A snic2018-8-151
#SBATCH -p core
Rscript --vanilla /proj/snic2018-8-151/private/Github/CIMseq.testing/inst/analysis/MGA.analysis_SI_strictClasses/scripts/CIMseqSwarm.R 36
