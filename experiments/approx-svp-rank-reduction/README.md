## Overview

This code implements the refined Approx-SVP rank reduction analysis. The project evaluates lattice security by analyzing the hardness of the Approximate Shortest Vector Problem (Approx-SVP) using advanced lattice reduction techniques and validates them through experiments on high-dimensional q-ary lattices.

## Structure

- **rank-reduction.py**: Main solver for q-ary lattice challenges
- **run_rank-reduction.sh**: Batch experiment runner
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
 chmod +x run_rank-reduction.sh
./run_rank-reduction.sh
```