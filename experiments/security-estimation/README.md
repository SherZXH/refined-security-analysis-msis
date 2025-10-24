# MSIS Security Estimation of Dilithium

## Overview

This code provides MSIS security analysis for the **Dilithium digital signature scheme**, a lattice-based post-quantum cryptographic algorithm selected as a finalist in the NIST post-quantum cryptography standardization process. The project evaluates security levels by analyzing the hardness of the Module Short Integer Solution (MSIS) problem using advanced lattice reduction techniques.

## Structure

- **Dilithium.py**: Main entry point for security analysis, defines parameter sets and runs evaluations
- **MSIS_security.py**: Core security evaluation engine implementing BKZ and sieve attack simulations  
- **model_BKZ.py**: BKZ lattice reduction algorithm modeling with classical, quantum, and plausible cost models
- **proba_util.py**: Probability calculation utilities for Gaussian and binomial distributions
- **README.md**: This documentation file

## Build and Run

### Prerequisites
- Python 3.6 or higher
- Standard math libraries (included with Python)

### Basic Execution
To run the security analysis with default parameters:
```bash
 python Dilithium.py
```

## Using "Dimension for Free"(D4f) Technique

###  About D4f Technique
Dimension-for-Free (D4f) is an advanced lattice reduction optimization technique proposed by Ducas [D4F2018]. The method involves initially sieving on a locally projected lattice $\mathcal{L}{[l, n]}$, then lifting all vectors in the list to $\mathcal{L}{[1, n]}$. If $\mathbf{v}$ represents a short vector in $\mathcal{L}{[1, n]}$, then $\mathbf{v}{[l, n]}$ is likely to also be short and thus contained in the sieving list. This allows sieving on $\mathcal{L}{[l, n]}$ to potentially recover the shortest vector in $\mathcal{L}{[1, n]}$, effectively providing $l$ dimensions "for free".

The asymptotic free dimension $l$ is given by $\frac{n \ln (4 / 3)}{\ln (n / 2 \pi e)}$, though practical implementations [D4F2018] can achieve up to $n / \ln (n)$.

###  To enable the D4f technique as used in the paper, uncomment the following lines in MSIS_security.py:

1. In the function SIS_linf_cost_our, around lines 196-197 (as per the provided code):
```bash
# bkz_cost, p_bkz = SIS_linf_cost_sub(i,j,l,B,q,round(b - theo_dim4free_fun_optimistic(b)),1)
# sieve_cost, p_sieve = SIS_linf_cost_sub(i,j,l,B,q,round(rank - theo_dim4free_fun_optimistic(rank)),1)
```

Replace them with:
```bash
bkz_cost, p_bkz = SIS_linf_cost_sub(i,j,l,B,q,round(b - theo_dim4free_fun_optimistic(b)),1)
sieve_cost, p_sieve = SIS_linf_cost_sub(i,j,l,B,q,round(rank - theo_dim4free_fun_optimistic(rank)),1)
```

2. Similarly, in the same function, around lines 236-237 (as per the provided code):
```bash
# bkz_cost_prime, p_bkz_prime = SIS_linf_cost_sub(i_prime,j_prime,l_prime,B,q,round(b_prime - theo_dim4free_fun_optimistic(b_prime)),1)
# sieve_cost_prime, p_sieve_prime = SIS_linf_cost_sub(i_prime,j_prime,l_prime,B,q,round(rank_prime - theo_dim4free_fun_optimistic(rank_prime)),1)
```

Replace them with:
```bash
bkz_cost_prime, p_bkz_prime = SIS_linf_cost_sub(i_prime,j_prime,l_prime,B,q,round(b_prime - theo_dim4free_fun_optimistic(b_prime)),1)
sieve_cost_prime, p_sieve_prime = SIS_linf_cost_sub(i_prime,j_prime,l_prime,B,q,round(rank_prime - theo_dim4free_fun_optimistic(rank_prime)),1)
```

3. To run the security analysis with D4f:
```bash
 python Dilithium.py
```
