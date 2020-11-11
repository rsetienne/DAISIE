#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --partition=gelifes
#SBATCH --mem=500000
Rscript test.R
