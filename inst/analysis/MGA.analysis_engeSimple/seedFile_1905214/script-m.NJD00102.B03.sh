#!/bin/bash -l
#SBATCH -J CIMseq-m.NJD00102.B03
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_engeSimple/seedFile_1905214_2/CIMseq-m.NJD00102.B03-%j.out
#SBATCH -t 00:30:00
#SBATCH -n 2
#SBATCH -A snic2018-8-151
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_engeSimple/scripts/CIMseqSwarm.R m.NJD00102.B03
