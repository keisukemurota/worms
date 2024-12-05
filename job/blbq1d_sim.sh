#!/bin/bash

# n: Calculate the total number of jobs for BLBQ model

calculate_total_jobs() {
    J0_values=(1)  # Define J0 values
    J1_values=($(seq -1 0.05 2))  # Define J1 values
    hz_values=(0)          # Define hz values
    hx_values=($(seq 0 0.05 1.05))          # Define hx values

    num_J0=${#J0_values[@]}
    num_J1=${#J1_values[@]}
    num_hz=${#hz_values[@]}
    num_hx=${#hx_values[@]}

    total_jobs=$((num_J0 * num_J1 * num_hz * num_hx))  # Updated total jobs calculation
}

# Function to calculate J0, J1, hx, and hz from a single integer input
calculate_parameters() {
    calculate_total_jobs
    local task_id=$1

    i=$((task_id / (num_J1 * num_hz * num_hx)))
    j=$((task_id % (num_J1 * num_hz * num_hx) / (num_hz * num_hx)))
    z=$((task_id % (num_hz * num_hx) / num_hx))
    x=$((task_id % num_hx))

    J0=$(printf "%.3f" "${J0_values[$i]}")
    J1=$(printf "%.3f" "${J1_values[$j]}")
    hz=$(printf "%.3f" "${hz_values[$z]}")
    hx=$(printf "%.3f" "${hx_values[$x]}")
}



echo_jobs() {

    # for loop for all parameter to confirm which parameter will be used
    for i in $(seq 0 $((total_jobs - 1)))
    do
        calculate_parameters "$i"
        echo "J0=${J0}, J1=${J1}, hz=${hz} and hx=${hx}"
    done
    echo "Total number of jobs for BLBQ model: $total_jobs"
}

# Core job function for BLBQ model
run_job() {
    local task_id=$1
    local project_dir=$2
    local n_cpu=$3

    calculate_parameters "$task_id"

    LT=1
    SWEEPS=1000000
    EPOCH=1000
    M=40
    model_name="BLBQ1D"
    log_dir="${project_dir}/job/log/${model_name}"
    [ ! -d "$log_dir" ] && mkdir -p "$log_dir" && echo "Created log directory $log_dir"
    log_file="${log_dir}/J0_${J0}_J1_${J1}_hz_${hz}_hx_${hx}_output.log"
    link_dir="${project_dir}/job/link/${model_name}"
    symbolic_link="${link_dir}/J0_${J0}_J1_${J1}_hz_${hz}_hx_${hx}_lt_${LT}"


    # Create the link directory if it does not exist
    [ ! -d "$link_dir" ] && mkdir -p "$link_dir" && echo "Created link directory $link_dir"
    echo "Symbolic link for BLBQ1D model: $symbolic_link"
    echo "Project directory: $project_dir"

    cd "$project_dir/python/rmsKit" || return

    # if n_cpu is positive, then use that number of CPUs
    if [ "$n_cpu" -gt 0 ]; then
        echo "Using $n_cpu CPUs"
        export OMP_NUM_THREADS=$n_cpu
        export MKL_NUM_THREADS=$n_cpu
    else
        echo "Using all available CPUs"
    fi
    source /opt/materiapps-intel/env.sh
    source ~/worms/myenv/bin/activate


    # Assuming the environment is already activated and required modules loaded
    python -u optimize_loc.py -m $model_name -o Adam -lr 0.001 -e $EPOCH -M $M -lt $LT \
        -J0 "$J0" -J1 "$J1" -hx "$hx" -hz "$hz" -n "$n_cpu" \
        --symoblic_link "$symbolic_link" \
        --stdout >> "$log_file"

    echo "Finished optimization for BLBQ model with J0=${J0}, J1=${J1}, hz=${hz} and hx=${hx} in CPU ${n_cpu}"


    python -u -m run_worm -m $model_name --path "$symbolic_link" -s $SWEEPS -n "$n_cpu" -k 10 --stdout  >> "$log_file"

    python -u -m run_worm -m $model_name --path "$symbolic_link" -s $SWEEPS --original -n "$n_cpu" --stdout  >> "$log_file"


    echo "Finished BLBQ model job with J0=${J0}, J1=${J1}, hz=${hz} and hx=${hx} in CPU ${n_cpu}"

    echo "Cleaning up existing symbolic links $symbolic_link"

    # Check if symbolic link already exists and remove it if it does
    # [ -L "$symbolic_link" ] && unlink "$symbolic_link"

}

# Initialize script
calculate_total_jobs
echo_jobs
echo "Sourcing script blbq-sim.sh"
export -f run_job
export -f calculate_parameters
export -f calculate_total_jobs
