#!/bin/bash -l
#SBATCH -J CIMseq-m.NJA01301.L24
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_allEngeDeconv/tmp/CIMseq-m.NJA01301.L24-%j.out
#SBATCH -t 72:00:00
#SBATCH -n 2
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_allEngeDeconv/scripts/CIMseqSwarm.R m.NJA01301.L24
