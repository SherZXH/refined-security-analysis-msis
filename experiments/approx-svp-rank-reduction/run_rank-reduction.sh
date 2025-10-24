#!/bin/bash

declare -A PARAM_SETS=(
      ["q850_exp_g1668906_beta136_practical"]="1.668906:../../data/approx-svp-rank-reduction/lattice-basis/pump_and_jumpbkz-blocksize-136-jump-003-tours-000-dim-850-cols-230.txt:164:230"
      ["q925_exp_g2777266_beta125_practical"]="2.777266:../../data/approx-svp-rank-reduction/lattice-basis/pump_and_jumpbkz-blocksize-125-jump-001-tours-000-dim-925-cols-295.txt:153:295"
      ["q1025_exp_g2891519_beta131_practical"]="2.891519:../../data/approx-svp-rank-reduction/lattice-basis/pump_and_jumpbkz-blocksize-140-jump-001-tours-000-dim-1025-cols-330.txt:170:330"  
      ["q850_exp_g1668906_beta136_compact"]="1.668906:../../data/approx-svp-rank-reduction/lattice-basis/pump_and_jumpbkz-blocksize-136-jump-003-tours-000-dim-850-cols-230.txt:149:230"
      ["q925_exp_g2777266_beta125_compact"]="2.777266:../../data/approx-svp-rank-reduction/lattice-basis/pump_and_jumpbkz-blocksize-125-jump-001-tours-000-dim-925-cols-295.txt:134:295"
      ["q1025_exp_g2891519_beta131_compact"]="2.891519:../../data/approx-svp-rank-reduction/lattice-basis/pump_and_jumpbkz-blocksize-140-jump-001-tours-000-dim-1025-cols-330.txt:162:330"
)

ORDERED_KEYS=(
    "q850_exp_g1668906_beta136_practical"
    "q925_exp_g2777266_beta125_practical"
    "q1025_exp_g2891519_beta131_practical"
    "q850_exp_g1668906_beta136_compact"
    "q925_exp_g2777266_beta125_compact"
    "q1025_exp_g2891519_beta131_compact"
)

NUMA_NODE_0_CPUS="0-35,72-107"
NUMA_NODE_1_CPUS="36-71,108-143"
THREADS_PER_NODE=72

# Create log directory (important: ensure directory exists)
mkdir -p logs

# Simplified GPU quick check (only run at start and end)
quick_gpu_check() {
    local exp_name=$1
    local stage=$2  # start or end
    
    {
        echo "=== GPU Quick Check ($stage) $(date) ==="
        nvidia-smi --query-gpu=index,name,memory.used,memory.total --format=csv
    } >> "logs/${exp_name}_quick_gpu.log"
}

for exp_name in "${ORDERED_KEYS[@]}"; do
    # Fix: add -r parameter to prevent backslash escaping
    IFS=':' read -r gamma filename rank cols <<< "${PARAM_SETS[$exp_name]}"
    
    if [[ ! -f $filename ]]; then
        echo "Error: File does not exist - $filename"
        echo "Problem occurred with parameter set: $exp_name"
        exit 1
    fi
    
    echo "================================"
    echo " Starting experiment set: $exp_name"
    echo " Gamma: $gamma | Preprocessed file: $filename"
    echo " Rank: $rank | Cols: $cols"
    echo "================================"
    
    for ((run=1; run<=1; run++)); do
        log_file="logs/${exp_name}_run.log"
        
        echo "$(date +'%Y-%m-%d %T') Starting run ${run}" | tee -a "$log_file"
        
        # Quick GPU check before starting
        quick_gpu_check "$exp_name" "start"
         
        # Execute main program
        echo "NUMA binding: CPU=$NUMA_NODE_1_CPUS, memory=node 1" | tee -a "$log_file"
        
        # Fix: Add error handling to ensure command execution
        if numactl --interleave=all python3 rank-reduction.py "$cols" \
            --gamma "$gamma" \
            --preprocessed-file "$filename" \
            --rank "$rank" \
            --threads "$THREADS_PER_NODE" \
            --gpus 4 2>&1 | tee -a "$log_file"
        then
            echo "Program executed successfully" | tee -a "$log_file"
        else
            echo "Program execution failed" | tee -a "$log_file"
        fi
        
        quick_gpu_check "$exp_name" "end"
        
        # Check exit status and provide suggestions
        exit_code=${PIPESTATUS[0]}
        if [ $exit_code -ne 0 ]; then
            echo "Warning: Program exited abnormally (code: $exit_code)" | tee -a "$log_file"
            echo "Suggested troubleshooting steps:" | tee -a "$log_file"
            echo "1. Check GPU memory usage" | tee -a "$log_file"
            echo "2. Check detailed error log: $log_file" | tee -a "$log_file"
            echo "3. Try reducing GPU count or thread count" | tee -a "$log_file"
            # Fix: Option to continue running next experiment instead of exiting
            # exit 1  # Uncomment this line if you want to stop the entire script on error
        fi
        
        echo "$(date +'%Y-%m-%d %T') Run ${run} completed" | tee -a "$log_file"
        echo "--------------------------------" >> "$log_file"
    done
done