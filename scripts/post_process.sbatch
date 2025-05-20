#!/bin/bash
#SBATCH --account=hai_solai
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=128
#SBATCH --job-name=jV_sim
#SBATCH --time=24:00:00
#SBATCH --partition=scraper
#SBATCH --nodes=1
source $HOME/.bashrc
echo "Running"

PROJECT_ROOT=/p/project1/solai/oestreicher1/repos/simsalabim/simulation_tool

cd ${PROJECT_ROOT}

source ${PROJECT_ROOT}/.venv/bin/activate

srun python ${PROJECT_ROOT}/scripts/post_process.py