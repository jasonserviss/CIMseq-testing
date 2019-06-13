#!/bin/bash -l
#SBATCH -J CIMseq-s.NJA02001.H06
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/deconvoluteSingletsMouse/seedFile_190610/CIMseq-s.NJA02001.H06-%j.out
#SBATCH -t 06:00:00
#SBATCH -n 2
#SBATCH -A snic2018-8-151
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/deconvoluteSingletsMouse/scripts/CIMseqSwarm.R s.NJA02001.H06 seedFile_190610
