#!/bin/bash

declare -A PARAM_SETS=(
    ["q1050_exp_g2892289_beta_116_pnjbkz_118"]="2.892289:../../data/lattice-basis-reduction/pump_and_jumpbkz-vs-wnjbkz/lattice-basis/pump_and_jumpbkz-116-jump-2-dim-1050-cols-340.txt:118:2:1:1:161:340"
    ["q1050_exp_g2892289_beta_118_pnjbkz_121"]="2.892289:../../data/lattice-basis-reduction/pump_and_jumpbkz-vs-wnjbkz/lattice-basis/pump_and_jumpbkz-118-jump-2-dim-1050-cols-340.txt:121:2:1:1:161:340"
    ["q1050_exp_g2892289_beta_118_wnjbkz_118"]="2.892289:../../data/lattice-basis-reduction/pump_and_jumpbkz-vs-wnjbkz/lattice-basis/pump_and_jumpbkz-118-jump-2-dim-1050-cols-340.txt:118:2:1:2:161:340"
    ["q1050_exp_g2892289_beta_118_wnjbkz_121"]="2.892289:../../data/lattice-basis-reduction/pump_and_jumpbkz-vs-wnjbkz/lattice-basis/pump_and_jumpbkz-118-jump-2-dim-1050-cols-340.txt:121:6:1:2:161:340"
    ["q1050_exp_g2892289_beta_109_headwnjbkz_111"]="2.892289:../../data/lattice-basis-reduction/wnjbkz-vs-headwnjbkz/lattice-basis/headwnjbkz-blocksize-109-jump-1-rank-161-dim-1050-cols-340.txt:111:1:1:3:161:340"
    ["q1050_exp_g2892289_beta_109_wnjbkz_111"]="2.892289:../../data/lattice-basis-reduction/wnjbkz-vs-headwnjbkz/lattice-basis/wnjbkz-blocksize-109-jump-1-dim-1050-cols-340.txt:111:1:1:2:161:340"
    ["q1050_exp_g2892289_beta_111_headwnjbkz_113"]="2.892289:../../data/lattice-basis-reduction/wnjbkz-vs-headwnjbkz/lattice-basis/headwnjbkz-blocksize-111-jump-1-rank-161-dim-1050-cols-340.txt:113:2:1:3:161:340"
    ["q1050_exp_g2892289_beta_111_wnjbkz_113"]="2.892289:../../data/lattice-basis-reduction/wnjbkz-vs-headwnjbkz/lattice-basis/wnjbkz-blocksize-111-jump-1-dim-1050-cols-340.txt:113:2:1:2:161:340"
    ["q1050_exp_g2892289_beta_113_headwnjbkz_118"]="2.892289:../../data/lattice-basis-reduction/wnjbkz-vs-headwnjbkz/lattice-basis/headwnjbkz-blocksize-113-jump-2-rank-161-dim-1050-cols-340.txt:118:2:1:3:161:340"
    ["q1050_exp_g2892289_beta_113_wnjbkz_118"]="2.892289:../../data/lattice-basis-reduction/wnjbkz-vs-headwnjbkz/lattice-basis/wnjbkz-blocksize-113-jump-2-dim-1050-cols-340.txt:118:2:1:2:161:340"
)

ORDERED_KEYS=(
    "q1050_exp_g2892289_beta_116_pnjbkz_118"
    "q1050_exp_g2892289_beta_118_pnjbkz_121"
    "q1050_exp_g2892289_beta_118_wnjbkz_118"
    "q1050_exp_g2892289_beta_118_wnjbkz_121"
    "q1050_exp_g2892289_beta_109_headwnjbkz_111"
    "q1050_exp_g2892289_beta_109_wnjbkz_111"
    "q1050_exp_g2892289_beta_111_headwnjbkz_113"
    "q1050_exp_g2892289_beta_111_wnjbkz_113"
    "q1050_exp_g2892289_beta_113_headwnjbkz_118"
    "q1050_exp_g2892289_beta_113_wnjbkz_118"
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
    IFS=':' read -r gamma filename blocksize jump tours mode rank cols <<< "${PARAM_SETS[$exp_name]}"
    
    if [[ ! -f $filename ]]; then
        echo "Error: File does not exist - $filename"
        echo "Problem occurred with parameter set: $exp_name"
        exit 1
    fi
    
    echo "================================"
    echo " Starting experiment set: $exp_name"
    echo " Gamma: $gamma | Preprocessed file: $filename"
    echo " Blocksize: $blocksize | jump: $jump |tours: $tours |mode: $mode"
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
        if numactl --interleave=all python3 wnjbkz.py "$cols" \
            --gamma "$gamma" \
            --preprocessed-file "$filename" \
            --blocksize "$blocksize" \
            --jump "$jump" \
            --tours "$tours" \
            --rank "$rank" \
            --mode "$mode" \
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