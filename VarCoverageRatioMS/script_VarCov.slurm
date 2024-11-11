#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -t 1-
#SBATCH --partition=day
#SBATCH -J var_coverage
#SBATCH -o ./log/log.var_coverage_%A_%a
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --array=1
#SBATCH --mem=64g

echo '-------------------------------'
cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo '-------------------------------'
echo -e '\n\n'

module load R

Rscript ./WR_Coverage.R $SLURM_ARRAY_TASK_ID
