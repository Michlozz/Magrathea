#!/bin/bash

## resource allocation
#SBATCH --job-name=planet_job
#SBATCH --output=planet_output.txt
#SBATCH --error=planet_error.txt 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1       # Number of tasks (CPU cores) per node
#SBATCH --cpus-per-task=1
#SBATCH --cpus-per-task=1         # Number of CPU threads per task
#SBATCH --mem=8g
#SBATCH --time=02:30:00           # Maximum execution time (hh:mm:ss)

#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --partition=work
#SBATCH --account=michaell

## modules and envs
module load  GCC  GSL

## runtime
./planet
