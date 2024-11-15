#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 7-
#SBATCH --partition=week
#SBATCH -J wr_quad_imb
#SBATCH -o ./log/log.wr_quad_imb_%A_%a
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

Rscript ./WR_Sim_quad_imb.R $SLURM_ARRAY_TASK_ID
