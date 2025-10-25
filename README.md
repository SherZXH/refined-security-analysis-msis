# Refined Security Analysis of MSIS via Improved Rank Reduction

## Overview

This repository contains the implementation and experimental code for the paper **"A Refined Security Analysis of MSIS via Improved Rank Reduction and Application on Dilithium"**. Our work presents improved Approx-SVP rank reduction conditions and novel lattice reduction algorithms (WnjBKZ and HeadWnjBKZ) that significantly advance the state-of-the-art in lattice-based cryptanalysis.

**Key Contributions:**
- Two refined Approx-SVP rank reduction conditions (practical and compact)
- WnjBKZ and HeadWnjBKZ algorithms for high-dimensional lattice reduction  
- Improved security estimation for MSIS problem in Dilithium
- Experimental validation showing up to 60× speedup on lattice challenges

## Environment Setup

### Prerequisite: GPU-G6K Installation

**Our algorithms are built upon the GPU-G6K library and require its proper installation first:**

```bash
# Clone the GPU-G6K repository
git clone https://github.com/WvanWoerden/G6K-GPU-Tensor
cd G6K-GPU-Tensor

# Follow the GPU-G6K installation instructions
# This typically involves:
# 1. Installing CUDA toolkit
# 2. Setting up GPU drivers
# 3. Building with GPU support
# 4. Installing Python bindings

# Refer to GPU-G6K's documentation for detailed installation steps
```

## Algorithm Integration

### After installing GPU-G6K, replace the following files in your g6k/algorithms/ directory with our improved versions:

```bash
# Copy our modified algorithm files to GPU-G6K installation
cp path/to/our/code/bkz.py /path/to/g6k/installation/g6k/algorithms/
cp path/to/our/code/pump.py /path/to/g6k/installation/g6k/algorithms/  
cp path/to/our/code/workout.py /path/to/g6k/installation/g6k/algorithms/
```
**Modified Algorithms Include:**
- Enhanced BKZ with WnjBKZ and HeadWnjBKZ algorithms
- WnjBKZ and HeadWnjBKZ algorithms for high-dimensional lattice reduction  
- Improved Pump procedure with optimized D4f integration
- Extended WorkOut algorithm with adaptive reduction scope



## Experimental Reproduction

### Chapter 3 & Section 5.1: Approx-SVP Rank Reduction

```bash
cd experiments/approx-svp-rank-reduction
# Follow README.md for detailed instructions
./run_rank-reduction.sh
```
- Validates practical and compact rank reduction conditions
- Tests dimensions 850, 925, 1025
- Compares with original Laarhoven-Mariano approach

### Chapter 4 & Section 5.2: Lattice Basis Reduction

```bash
cd experiments/lattice-basis-reduction
# Follow README.md for detailed instructions  
./run_wnjbkz.sh
```
- WnjBKZ and HeadWnjBKZ algorithms for high-dimensional lattice reduction  
- Tests various blocksizes and jump parameters
- Evaluates saturation ratios and basis quality

### Chapter 6: Security Estimation

```bash
cd experiments/security-estimation
# Follow README.md for detailed instructions
python security_analysis.py
```
**Modified Algorithms Include:**
- Implements improved MSIS security model
- Analyzes Dilithium parameter sets  
- Compares with Core-SVP baseline

## Contributing
-For questions about code implementation or experimental reproduction:

-Check the detailed README files in each experiment subdirectory

-Review the logs in logs/ for expected output formats

-Examine data/ for input/output examples

-Contact the authors for specific implementation questions


