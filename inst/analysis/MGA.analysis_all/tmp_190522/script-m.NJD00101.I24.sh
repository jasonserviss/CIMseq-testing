#!/bin/bash -l
#SBATCH -J CIMseq-m.NJD00101.I24
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_all/tmp/CIMseq-m.NJD00101.I24-%j.out
#SBATCH -t 06:00:00
#SBATCH -n 2
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_all/scripts/CIMseqSwarm.R m.NJD00101.I24
