#!/bin/bash
#SBATCH --job-name run_mad_yg_train_1
#SBATCH --time 12:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 52
#SBATCH -o output_mad_yg_train_1.%j
#SBATCH -e error_mad_yg_train_1.%j

module load anaconda
conda activate planet_processing_test
python 8_run_mad_calibration_in_subprocess.py

