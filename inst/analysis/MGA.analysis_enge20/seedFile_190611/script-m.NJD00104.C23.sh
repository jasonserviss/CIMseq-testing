#!/bin/bash -l
#SBATCH -J CIMseq-m.NJD00104.C23
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_enge20/seedFile_190611/CIMseq-m.NJD00104.C23-%j.out
#SBATCH -t 48:00:00
#SBATCH -n 3
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_enge20/scripts/CIMseqSwarm.R m.NJD00104.C23 seedFile_190611
