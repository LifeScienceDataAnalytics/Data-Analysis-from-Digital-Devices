#!/bin/bash
#SBATCH -p draco2
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --array=0-11
#SBATCH --time=7-00:00:00
#SBATCH --output=/home/msood/ADNIVAMBN_paper_all/GridSearch/slurm.%A_%a.out
#SBATCH --error=/home/msood/ADNIVAMBN_paper_all/GridSearch/slurm.%A_%a.err


module load Anaconda3 foss/2020a
# module --ignore-cache load CUDA/10.1.168
source activate /home/msood/.conda/envs/tensorflow

python3 -u GridSearch_Altoida.py $SLURM_ARRAY_TASK_ID
