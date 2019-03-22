#!/bin/bash -l
#SBATCH -J CIMseq-m.NJD00104.E17
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_colonProximal/tmp/CIMseq-m.NJD00104.E17-%j.out
#SBATCH -t 72:00:00
#SBATCH -n 1
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_colonProximal/scripts/CIMseqSwarm.R m.NJD00104.E17
