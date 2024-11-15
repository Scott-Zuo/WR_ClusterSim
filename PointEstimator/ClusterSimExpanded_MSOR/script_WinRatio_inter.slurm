#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 7-
#SBATCH --partition=week
#SBATCH -J wr_interact_msor
#SBATCH -o ./log/log.wr_interact_msor_%A_%a
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

Rscript ./Sim_inter_parall_MSOR.R $SLURM_ARRAY_TASK_ID
