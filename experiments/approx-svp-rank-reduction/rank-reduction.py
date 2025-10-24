#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Q-ary Lattice Challenge Solver Command Line Client
Solves q-ary lattice challenges by finding vectors with norm less than given threshold
"""

import sys
import os

# Add the path to g6k directory (two levels up from current file)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

import time
import copy
import logging
import pickle as pickler
from collections import OrderedDict
from fpylll.util import gaussian_heuristic
from g6k.algorithms.workout import workout
from g6k.siever import Siever
from g6k.utils.cli import parse_args, run_all, pop_prefixed_params
from g6k.utils.stats import SieveTreeTracer
from g6k.utils.util import load_svpchallenge_and_randomize, load_matrix_file, db_stats
from fpylll.tools.bkz_stats import dummy_tracer
from math import log, exp, lgamma, pi, sqrt, e

from fpylll.algorithms.bkz2 import BKZReduction
from fpylll.tools.quality import basis_quality
from numpy import True_
from fpylll import BKZ as fplll_bkz, GSO, IntegerMatrix, LLL

from g6k.algorithms.bkz import pump_n_jump_bkz_tour, workout_bkz_tour, theo_dim4free_fun2
from g6k.algorithms.pump import pump



import numpy as np
from fpylll import FPLLL


def qary_kernel(arg0, params=None, seed=None):
    """
    Kernel function for solving q-ary lattice challenges
    Finds vectors with norm less than threshold defined by gamma parameter
    """
    logger = logging.getLogger('qary')

    # Handle single parameter case for Pool.map
    if params is None and seed is None:
        n, params, seed = arg0
    else:
        n = arg0

    params = copy.copy(params)

    # Extract various parameters
    load_matrix = params.pop("load_matrix")
    pre_bkz = params.pop("pre_bkz")
    pump_params = pop_prefixed_params("pump", params)
    workout_params = pop_prefixed_params("workout", params)
    verbose = params.pop("verbose")
    if verbose:
        workout_params["verbose"] = True
    challenge_seed = params.pop("challenge_seed")
    high_prec = params.pop("high_prec")
    trace = params.pop("trace")
    
    # Extract q-ary lattice specific parameters
    filename = params.pop("preprocessed_file")
    gamma = params.pop("gamma")
    rank = params.pop("rank")

    # Check if preprocessed file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Preprocessed matrix file not found: {filename}")

    # Load lattice basis from file
    A = IntegerMatrix.from_file(filename)
    g6k = Siever(A, params, seed=seed)

    print("GSO precision: ", g6k.M.float_type)

    # Initialize tracer for performance monitoring
    if trace:
        tracer = SieveTreeTracer(g6k, root_label=("qary-challenge", n), start_clocks=True)
    else:
        tracer = dummy_tracer
    
    # Apply LLL reduction to the basis
    g6k.lll(0, g6k.full_n)
    
    # Calculate Gaussian heuristic for the lattice
    gh_2 = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])

    # Calculate volume and related metrics
    r = [g6k.M.get_r(i, i) for i in range(n)]
    n = len(list(r))
    log_vol = sum([log(x) for x in r])
    log_vol_n = 1./n * (log_vol)
    vol_n = exp(log_vol_n)
    print("vol_n:%f" % sqrt(vol_n))

    # Calculate target norm based on gamma parameter
    goal_r0 = (gamma**2) * gh_2  # Use gamma from command line
    
    if verbose:
        print("gh = %f, goal_r0/gh = %f, r0/gh = %f, rf = %f" % (
            sqrt(gh_2), 
            sqrt(goal_r0/gh_2), 
            sqrt(sum([x*x for x in A[0]])/gh_2),
            sqrt(sum([x*x for x in A[0]])/vol_n)
        ))
        slope = basis_quality(g6k.M)["/"]
        print("Current Slope = %.5f\n" % slope)
        print("goal:%f" % (sqrt(goal_r0)))

    # Attempt to find short vector using pump algorithm
    for cnt in range(10):
        flast = pump(g6k, tracer, 0, rank, round(rank/log(rank)), goal_r0=goal_r0, verbose=True)[0]
        print("---------current_initial_vector:%.5f--------" % (g6k.M.get_r(0, 0)**0.5))
        
        # Check if we found a vector with norm less than target
        if g6k.M.get_r(0, 0)**0.5 < sqrt(goal_r0):
            # Ensure output directory exists before saving
            output_dir = "qarychallenge"
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                print(f"Created directory: {output_dir}")
            
            # Save successful result
            block_filename = "qarychallenge/pump-%02d-dim-%03d.txt" % (rank, g6k.full_n)
            fn = open(block_filename, "w")
            fn.write(str(g6k.M.B))
            fn.close()
            print(f"Saved result to: {block_filename}")
            break
    
    if verbose:
        logger.info("sol %d, %s" % (n, A[0]))

    # Calculate and log final norm information
    norm = sum([x*x for x in A[0]])
    if verbose:
        logger.info("norm %.1f ,hf %.5f" % (norm**.5, (norm/gh_2)**.5))

    tracer.exit()

    # Return statistics if tracing is enabled
    if hasattr(tracer, "trace"):
        stat = tracer.trace
        stat.data["flast"] = flast
        return stat
    else:
        return None


def qary():
    """
    Main function for solving q-ary lattice challenges
    Runs multiple experiments with different parameters
    """
    description = qary.__doc__

    # Parse command line arguments with default values
    args, all_params = parse_args(description,
                                  load_matrix=None,
                                  pre_bkz=None,
                                  verbose=True,
                                  challenge_seed=0,
                                  workout__dim4free_dec=2,
                                  trace=True,
                                  gamma=1.40,          # Default gamma value
                                  rank=100,
                                  preprocessed_file="pnjbkz-beta-136-jump-003-tours-000-dim-850-cols-230.txt")

    # Run experiments with all parameter combinations
    stats = run_all(qary_kernel, all_params.values(),
                    lower_bound=args.lower_bound,
                    upper_bound=args.upper_bound,
                    step_size=args.step_size,
                    trials=args.trials,
                    workers=args.workers,
                    seed=args.seed)

    # Create inverse mapping for parameter identification
    inverse_all_params = OrderedDict([(v, k) for (k, v) in all_params.items()])

    # Process and display results
    for (n, params) in stats:
        stat = stats[(n, params)]
        if stat[0] is None:
            logging.info("Trace disabled")
            continue

        if len(stat) > 0:
            # Calculate average statistics across trials
            cputime = sum([float(node["cputime"]) for node in stat])/len(stat)
            walltime = sum([float(node["walltime"]) for node in stat])/len(stat)
            flast = sum([float(node["flast"]) for node in stat])/len(stat)
            avr_db, max_db = db_stats(stat)
            
            # Format and print results
            fmt = "%48s :: m: %1d, n: %2d, cputime :%7.4fs, walltime :%7.4fs, flast : %2.2f, avr_max db: 2^%2.2f, max_max db: 2^%2.2f"
            print(fmt % (inverse_all_params[params], params.threads, n, cputime, walltime, flast, avr_db, max_db))
        else:
            logging.info("Trace disabled")

    # Save results to pickle file if requested
    if args.pickle:
        pickler.dump(stats, open("hkz-qary-%d-%d-%d-%d.sobj" %
                                 (args.lower_bound, args.upper_bound, args.step_size, args.trials), "wb"))


if __name__ == '__main__':
    qary()