#!/bin/bash -l
#SBATCH -J CIMseq-m.NJA01202.B23
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_SI/seedFile_190613/CIMseq-m.NJA01202.B23-%j.out
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_SI/scripts/CIMseqSwarm.R m.NJA01202.B23 seedFile_190613
