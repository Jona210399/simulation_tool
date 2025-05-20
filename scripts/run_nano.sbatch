#!/bin/bash
#SBATCH --cpus-per-task=128
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --tasks-per-node=1
#SBATCH --job-name=jV_sim
#SBATCH --time=1:00:00
#SBATCH --partition=short

source ~/.bashrc
echo "Running simulation..."

PROJECT_ROOT=/home/ws/gv6569/repos/simsalabim/simulation_tool

cd ${PROJECT_ROOT}

source ${PROJECT_ROOT}/.venv/bin/activate

python ${PROJECT_ROOT}/scripts/run.py