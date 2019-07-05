#!/bin/bash -l
#SBATCH -J CIMseq-m.NJA01202.P21
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_SI20/seedFile_190625/CIMseq-m.NJA01202.P21-%j.out
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_SI20/scripts/CIMseqSwarm.R m.NJA01202.P21 seedFile_190625
