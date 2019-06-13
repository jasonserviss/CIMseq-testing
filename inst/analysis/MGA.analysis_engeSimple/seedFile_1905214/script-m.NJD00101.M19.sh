#!/bin/bash -l
#SBATCH -J CIMseq-m.NJD00101.M19
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_engeSimple/tmp/CIMseq-m.NJD00101.M19-%j.out
#SBATCH -t 00:20:00
#SBATCH -n 2
#SBATCH -A snic2018-8-151
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_engeSimple/scripts/CIMseqSwarm.R m.NJD00101.M19
