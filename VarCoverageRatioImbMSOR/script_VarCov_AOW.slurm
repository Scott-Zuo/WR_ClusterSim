#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -t 3-
#SBATCH --partition=week
#SBATCH -J var_coverage_AOW_ImbMSOR
#SBATCH -o ./log/log.var_coverage_AOW_ImbMSOR_%A_%a
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

Rscript ./WR_Coverage_AOW.R $SLURM_ARRAY_TASK_ID
