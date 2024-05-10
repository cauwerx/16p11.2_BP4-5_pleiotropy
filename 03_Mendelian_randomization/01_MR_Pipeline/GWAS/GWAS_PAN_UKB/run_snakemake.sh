#!/bin/bash

#SBATCH --job-name='sn_GWAS'
#SBATCH --nodes=1
#SBATCH --time=3:00:00
#SBATCH --partition urblauna
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=1GB

snakemake -j 3 --cluster "sbatch -p urblauna --job-name='GWASprep' --time=2:00:00 --nodes=1 --cpus-per-task=1 --mem-per-cpu=40G"

