#!/bin/bash
#SBATCH --job-name="mers_sim"
#SBATCH --output="mers.%j.%N.out"
#SBATCH --partition=gpuA100x4 # gpu type
#SBATCH --mem=64G
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8   # spread out to use 1 core per numa, set to 64 if tasks is 1
#SBATCH --constraint="scratch"
#SBATCH --gpus-per-node=4 # number of gpus per node
#SBATCH --gpu-bind=closest   # select a cpu close to gpu on pci bus topology
#SBATCH --account=<account-name>    # match to a "Project" returned by the "accounts" command
#SBATCH --exclusive  # dedicated node for this job
#SBATCH --no-requeue
#SBATCH -t 18:00:00 # job timeout duration
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

MSBS_ROOT=/u/pkumar9/xlab/msbs

# Activate Mamba-forge.
eval "$(${MSBS_ROOT}/mambaforge/bin/conda shell.bash hook)"
source ${MSBS_ROOT}/mambaforge/etc/profile.d/mamba.sh

mamba activate openmm

python main.py 4l72 2 1
