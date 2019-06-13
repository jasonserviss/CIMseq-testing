#!/bin/bash -l
#SBATCH -J CIMseq-m.NJA01205.P06
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_allEngeDeconv/seedFile_190524/CIMseq-m.NJA01205.P06-%j.out
#SBATCH -t 06:00:00
#SBATCH -n 1
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_allEngeDeconv/scripts/CIMseqSwarm.R m.NJA01205.P06 seedFile_190524
