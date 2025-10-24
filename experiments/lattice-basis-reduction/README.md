# WnjBKZ and HeadWnjBKZ Implementation

## Overview

This code implements the **WnjBKZ** (WorkOut and jump BKZ) and **HeadWnjBKZ** algorithms for advanced lattice basis reduction. The project evaluates and validates these improved lattice reduction techniques through comprehensive experiments on high-dimensional q-ary lattices, addressing critical limitations in traditional approaches.

## Structure

- **wnjbkz.py**: Main entry point for WnjBKZ and HeadWnjBKZ algorithms
- **run_wnjbkz.sh**: Batch experiment runner for algorithm comparison  
- **README.md**: This documentation file

## Build and Run

### Prerequisites
- Python 3.7 or higher
- [GPU-G6K](https://github.com/WvanWoerden/G6K-GPU-Tensor) lattice reduction library (GPU-accelerated version)
- [FPLLL](https://github.com/fplll/fplll) library
- CUDA toolkit (compatible with your GPU)
- NumPy
- NVIDIA GPUs (required for GPU acceleration)

### Batch Experiments
To run all experiments from the paper:
```bash
chmod +x run_wnjbkz.sh
./run_wnjbkz.sh