#!/bin/bash -l
#SBATCH -J CIMseq-m.NJD00102.J03
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_enge/seedFile_190611/CIMseq-m.NJD00102.J03-%j.out
#SBATCH -t 06:00:00
#SBATCH -n 2
#SBATCH -A snic2018-8-151
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_enge/scripts/CIMseqSwarm.R m.NJD00102.J03 seedFile_190611