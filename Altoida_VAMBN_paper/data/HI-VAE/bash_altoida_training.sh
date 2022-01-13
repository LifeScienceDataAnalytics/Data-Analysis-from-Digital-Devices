#!/bin/bash
#SBATCH -J AltTraining
#SBATCH -n 192
#SBATCH -t 48:00:00
#SBATCH -p smp





# SBATCH -p 
# SBATCH -n 1
# SBATCH --cpus-per-task=16
# SBATCH --array=0-20
# SBATCH --time=7-00:00:00
# SBATCH --output=/home/msood/ParkVAMBNFinal/GridSearch/slurm.%A_%a.out
# SBATCH --error=/home/msood/ParkVAMBNFinal/GridSearch/slurm.%A_%a.err

module load Anaconda3 foss/2020a
# module --ignore-cache load CUDA/10.1.168
source activate /home/msood/.conda/envs/tensorflow

python3 -u 1_Altoida_HIVAE_training.py


