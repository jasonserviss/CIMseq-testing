#!/bin/bash -l
#SBATCH -J CIMseq-m.NJA01501.L23
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_manualColonSplit/tmp/CIMseq-m.NJA01501.L23-%j.out
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_manualColonSplit/scripts/CIMseqSwarm.R m.NJA01501.L23