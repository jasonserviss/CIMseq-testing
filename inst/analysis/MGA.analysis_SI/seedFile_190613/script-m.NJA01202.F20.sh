#!/bin/bash -l
#SBATCH -J CIMseq-m.NJA01202.F20
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_SI/seedFile_190613_2/CIMseq-m.NJA01202.F20-%j.out
#SBATCH -t 12:00:00
#SBATCH -n 2
#SBATCH -A snic2018-8-151
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_SI/scripts/CIMseqSwarm.R m.NJA01202.F20 seedFile_190613_2
