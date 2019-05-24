#!/bin/bash -l
#SBATCH -J CIMseq-m.NJD00103.J12
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_engeSimple/tmp/CIMseq-m.NJD00103.J12-%j.out
#SBATCH -t 00:30:00
#SBATCH -n 2
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_engeSimple/scripts/CIMseqSwarm.R m.NJD00103.J12
