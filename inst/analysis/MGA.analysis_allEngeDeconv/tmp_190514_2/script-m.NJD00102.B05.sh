#!/bin/bash -l
#SBATCH -J CIMseq-m.NJD00102.B05
#SBATCH -o /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_allEngeDeconv/tmp/CIMseq-m.NJD00102.B05-%j.out
#SBATCH -t 6:00:00
#SBATCH -n 2
#SBATCH -A snic2019-3-84
#SBATCH -p core
Rscript --vanilla /home/jason/Github/CIMseq.testing/inst/analysis/MGA.analysis_allEngeDeconv/scripts/CIMseqSwarm.R m.NJD00102.B05
