#!/bin/bash -l
#SBATCH -J CIMseq-m.NJA01401.A20
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_allEngeDeconv/seedFile_190524/CIMseq-m.NJA01401.A20-%j.out
#SBATCH -t 06:00:00
#SBATCH -n 1
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_allEngeDeconv/scripts/CIMseqSwarm.R m.NJA01401.A20 seedFile_190524
