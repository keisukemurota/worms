#!/bin/bash

#SBATCH --job-name=blbq1d_sim
#SBATCH --output=job_log/logfile_%A_%a.log
#SBATCH --array=0-128  # Adjust this for the number of jobs and parallelism
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --partition=cpu

# Source the blbl.sh script
source blbq1d_sim.sh
#source blbq1d_dimer_sim.sh
# source kh2d_sim.sh
#source hxyz1d_sim.sh
#source hxyz2d_sim.sh
# source ss2d_sim.sh
# source mg1d_sim.sh

# Set the project directory
PROJECT_DIR=$(dirname "$(pwd)")

# This environment variable is automatically set by Slurm to the current array index
SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID

n_job_per_run=11  # Adjust this for the number of jobs per run
# n_job_per_run=4  # Adjust this for the number of jobs per run

# Run the job (for loop)
# rm "$PROJECT_DIR"/job_log/log/*

for i in $(seq 0 $((n_job_per_run - 1))); do
    job_id=$((SLURM_ARRAY_TASK_ID * n_job_per_run + i)) 
    # n: if job_id exceeds the total number of jobs, then break
    if [ "$job_id" -ge "$total_jobs" ]; then
        break
        echo "Job $job_id exceeds the total number of jobs $total_jobs"
    fi
    echo "Running job $job_id"
    run_job "$job_id" "$PROJECT_DIR" -1
done
