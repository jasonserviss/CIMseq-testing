#!/bin/bash -l
#SBATCH -J CIMseq-m.NJA01303.F22
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_enge/seedFile_190628/CIMseq-m.NJA01303.F22-%j.out
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_enge/scripts/CIMseqSwarm.R m.NJA01303.F22 seedFile_190628