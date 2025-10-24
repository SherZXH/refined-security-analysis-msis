#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SVP Challenge Solver Command Line Client
"""
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
from math import log, exp, lgamma, pi,sqrt,e
import os

from fpylll.algorithms.bkz2 import BKZReduction
from fpylll.tools.quality import basis_quality
from numpy import True_
from fpylll import BKZ as fplll_bkz, GSO, IntegerMatrix, LLL

from g6k.algorithms.bkz import pump_n_jump_bkz_tour,workout_bkz_tour,theo_dim4free_fun2
from g6k.algorithms.pump import pump

import numpy as np
from fpylll import FPLLL

def dim4free_set(dimension,key,slope,rhf):
    """
    Return number of dimensions for free according to different situation.

    :param dimension: the dimension of pump sieving algorithm
    :param key: 
        0-> (exact-SVP) default g6k, expected number of dimensions for free from exact-SVP experiments
        1-> (exact-SVP) ducas18, optimistic number of dimensions for free from hkz-reduced (bkz-n-reduced)
    """
    if key == 0:
        return int(11.5 + 0.075*dimension)
    elif key == 1:
        return int(dimension*log(4/3.)/log(dimension/2./pi/e))
    elif key == 2:
        return int(dimension*log(4/3.)/log(dimension/2./pi))
    elif key == 3:
        return int(log(sqrt(4/3.))/log(rhf))
    elif key == 4:
        return int(log(sqrt(4/3.))*(-2*dimension)/((dimension-1)*slope))
    elif key == 5:
        for f in range(1,int(dimension/2)):
            right = log(sqrt((4*(dimension-f))/(3*dimension)))/log(rhf)
            if f >=right:
                return f
    elif key == 6:
        for f in range(1,int(dimension/2)):
            right = log(sqrt(4/3))/log(rhf)
            if f >=right:
                return f
    elif key == 7:
        for f in range(1,int(dimension/2)):
            right = log(sqrt((4*(dimension-f))/(3*dimension)))*(-2*dimension)/((dimension-1)*slope)
            if f >=right:
                return f



def asvp_kernel(arg0,params=None, seed=None):
    logger = logging.getLogger('asvp')

    # Pool.map only supports a single parameter
    if params is None and seed is None:
        n, params, seed = arg0
    else:
        n = arg0

    params = copy.copy(params)

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
    
    """
    if load_matrix is None:
        A, bkz = load_svpchallenge_and_randomize(n, s=challenge_seed, seed=seed)
        if verbose:
            print("Loaded challenge dim %d" % n)
        if pre_bkz is not None:
            par = BKZ_FPYLLL.Param(pre_bkz, strategies=BKZ_FPYLLL.DEFAULT_STRATEGY, max_loops=1)
            bkz(par)

    else:
        A, _ = load_matrix_file(load_matrix, doLLL=False, high_prec=high_prec)
        if verbose:
            print("Loaded file '%s'" % load_matrix)
    """
    # 替换原来的filename和goal_r0计算部分
    filename = params.pop("preprocessed_file")
    gamma = params.pop("gamma")


    # 自动校正新旧命名差异
    filename = (
    filename
    .replace("pump_and_jumpbkz", "pnjbkz")  # 将长名称转为实际使用的缩写
    .replace("Tpump_and_jumpbkz", "Tpnjbkz")
    .replace("blocksize", "beta")
    )

    print(f"\n🔧 自动校正后的文件名: {filename}")

    # 检查文件存在
    if not os.path.exists(filename):
        raise FileNotFoundError(f"预处理矩阵文件不存在: {filename}")


    #filename ="svpchallenge/120-midmat-bkz49-pnjbkz-beta65-jump1-reduced-hf-1.380618-T(bkz)-20.071112s-T(pnjbkz)-1025.015023s.txt"
    #A,_= load_matrix_file(filename,randomize=False,high_prec=True)

    A = IntegerMatrix.from_file(filename)
    g6k = Siever(A, params, seed=seed)

    print("GSO precision: ", g6k.M.float_type)

    if trace:
        tracer = SieveTreeTracer(g6k, root_label=("svp-challenge", n), start_clocks=True)
    else:
        tracer = dummy_tracer
    g6k.lll(0, g6k.full_n)
    gh_2 = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])

    r = [g6k.M.get_r(i, i) for i in range(n)]
    n = len(list(r))
    log_vol = sum([log(x) for x in r])
    log_vol_n =  1./n * (log_vol)
    vol_n = exp(log_vol_n)
    print("vol_n:%f"%sqrt(vol_n))
   
    """
    gamma = params.pop("gamma")
   
    if gamma is None:
        goal_r0 = (1.14**2) * gh_2
    else:
        goal_r0 = gamma**2 * gh_2
    """
    # 修改goal_r0计算方式
    goal_r0 = (gamma**2) *   gh_2  # 使用从命令行传入的gamma值# 使用从命令行传入的gamma值
    if verbose:
        print("gh = %f, goal_r0/gh = %f, r0/gh = %f, rf = %f" % (sqrt(gh_2), sqrt(goal_r0/gh_2), sqrt(sum([x*x for x in A[0]])/gh_2),sqrt(sum([x*x for x in A[0]])/vol_n)))
        slope = basis_quality(g6k.M)["/"]
        print("Current Slope = %.5f\n" % slope)
    
    """
    prebeta = 59

    for blocksize in range(10,prebeta+1):
        print("Preprocess: Starting a BKZ-%d tour. " % (blocksize))
        bkz = BKZReduction(g6k.M)
        par = fplll_bkz.Param(blocksize,strategies=fplll_bkz.DEFAULT_STRATEGY,max_loops=1)
        bkz(par)
        block_filename = "qarychallenge/dim-%d/prebkz-%02d-dim-%03d.txt" % (n,blocksize, n)
        fn = open(block_filename, "w")
        fn.write(str(A))
        fn.close()
        slope = basis_quality(g6k.M)["/"]
        print("Current Slope = %.5f\n" % slope)
    """


    
    #beta_list =  [(121, 2, 2), (123, 2, 1), (125, 2, 2), (128, 2, 2), (130, 2, 1), (133, 2, 2), (135, 2, 1), (137, 2, 1)]
    
    #beta_list =  [(127, 2,1),(131, 2,1)]
    
    '''
    beta_list =   [(132, 7, 2), (136, 7, 1)]
    beta_list = [list(item) for item in  beta_list]

    print(beta_list)
    print(pump_params)
    
    for seq in beta_list:
            T0 = time.time()
            for t in range(seq[2]):
                T1 =time.time() 
                print("Starting a %dth workoutBKZ-%d-%d tour." % (t,seq[0],seq[1]))
                #pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize=beta,  jump = jump, verbose=True, extra_dim4free=0, pump_params=pump_params,goal_r0=goal_r0)
                workout_bkz_tour(g6k, dummy_tracer, blocksize=seq[0], extra_dim4free=0,jump = seq[1],tours = t,total_tours = seq[2],pump_params=pump_params,goal_r0=goal_r0) 
                print("tourtime: %f sec." %(time.time() - T1)) 
                print("||b0|| = %f. \n" % g6k.M.get_r(0, 0)**0.5)
                slope = basis_quality(g6k.M)["/"]
                print("Current Slope = %.5f\n" % slope)
                block_filename = "qarychallenge/dim-%d/wnjbkz-%d-jump-%d-%dth-tours.txt" % (n,seq[0], seq[1],t)
                #block_filename = "svpchallenge/dim-%d/wnjbkz-%d-jump-%d-%dth-tours.txt" % (n,seq[0], seq[1],t+1)
                fn = open(block_filename, "w")
                fn.write(str(g6k.M.B))
                fn.close()
            T_BKZ = time.time() - T0
            print("walltime: %f sec." % T_BKZ)
    
    g6k.lll(0, g6k.full_n)
    print("After workout-BKZ: norm %.5f" % (g6k.M.get_r(0, 0)**0.5))
    '''
    
    '''
    rank = 177
    flast = pump(g6k, tracer, 0,rank, round(rank/log(rank)),goal_r0=goal_r0, verbose = True)[0]
    block_filename = "svpchallenge/pump-%02d-dim-%03d.txt" % (rank,g6k.full_n)
    fn = open(block_filename, "w")
    fn.write(str(g6k.M.B))
    fn.close()
    '''


    '''
    for sample_dimension in range(80,180,1):
        #gh = (4.0/3.0)*gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n-sample_dimension,n)])
        gh = (4.0/3.0)*gaussian_heuristic([g6k.M.get_r(i, i) for i in range(sample_dimension)])
        if(gh<goal_r0):
            break
    '''
   
    
    #######rank reduction#####
    
    q = g6k.M.get_r(0, 0)
    print("initial_vector:%f"%(g6k.M.get_r(0, 0) ** 0.5))

    r = [g6k.M.get_r(j, j) for j in range(n)]
    log_vol = sum([log(x) for x in r])
    log_vol_n = 1.0 / n * log_vol
    vol_n = sqrt(exp(log_vol_n))
    gh = gaussian_heuristic(r) ** 0.5
    
    # 计算初始向量长度和RHF
    hf = g6k.M.get_r(0,0)**0.5/vol_n
    initial_vector = g6k.M.get_r(0, 0) ** 0.5
    initial_rhf = pow(hf,1.0/n)
    
    sample_dimension_list = [160,156,164]
    gamma_list = [2.8915,2.8565,2.8390]



    print(f"sample_dimension: %d"%(sample_dimension_list[0]))

    current_slope = basis_quality(g6k.M)["/"]
    print(f"initial_basis_slope: {current_slope/2:.5f}")
    for rank in sample_dimension_list:
        gamma = gamma_list[int((rank - sample_dimension_list[0])/4)]
        goal_r0 = (gamma**2) * gh_2
        gh_rank = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(rank)])

        print("Gaussian heuristic: %f" % sqrt(gh_rank))
        sigma = (((2*sqrt(3.0)+3)/6) * sqrt(gh_rank))/sqrt(g6k.full_n)
        print("standard deviation: %f"%sigma)
    

        T0 = time.time()
        max_goal = (4.0/3.0)*gaussian_heuristic([g6k.M.get_r(i, i) for i in range(rank)])
        print("goal:%f"%(sqrt(goal_r0)))
        print("max_goal:%f"%(sqrt(max_goal)))

        #dim4free_max = dim4free_set(rank,6,current_slope,initial_rhf)


        #print("dim4free:%d"%dim4free_max)

        #sample 3.2*(2^(0.2075*#sample_dimension)) vectors and select 1000 vectors whose length is shorter than 0.99*||b1||
        #pump(g6k, dummy_tracer,0,n,n-sample_dimension, verbose=True, down_sieve = True,goal_r0=goal) 
        flast = pump(g6k, tracer, 0,rank,round(rank/log(rank)),goal_r0=goal_r0, down_sieve = True,verbose = True)[0]
        #flast = pump(g6k, tracer, 0,rank,dim4free_max,goal_r0=goal_r0, down_sieve = True,verbose = True)[0]
        

        
        print("min_current_norm:%f"%(sqrt(g6k.M.get_r(0, 0))))
        #######rank reduction#####
    
        
    '''
    T_SAMPLE = time.time() - T0
    print("walltime: %f sec." %T_SAMPLE)
    print("sample finished.")

    T1 = time.time()
    db = list(g6k.itervalues())
    print("the size of database is %d"%(len(db)))
    found = 0

    sample_vector_list= []
    vector_norm_distribution = []
    min_vector_norm = float("inf")

    not_found = 0

    v = A.multiply_left(db[1000])
    print(v)
    l = sum(v_**2 for v_ in v)
    print(l)

    w = g6k.M.from_canonical(v)
    print(w)

    l_c = sum(w_**2*g6k.M.get_r(i, i) for i,w_ in enumerate(w))
    print(l_c)

    w_c = [w_*(g6k.M.get_r(i, i)**0.5) for i,w_ in enumerate(w)]
    print(w_c)

    l_w_c = sum(v_**2 for v_ in w_c)
    print(l_w_c)



    

    for j in range(int(min(1e6,len(db)))):
    #for j in range(len(db)):
        v = A.multiply_left(db[j])
        l = sum(v_**2 for v_ in v)

        if l <  1.0*max_goal and l!=q:
            #print(sqrt(l/max_goal) , v)
            found += 1
            not_found = 0
            sample_vector_list.append(sqrt(l))
            if l < min_vector_norm:
                min_vector_norm = l
            w = g6k.M.from_canonical(v)
            w_c = [w_*(g6k.M.get_r(i, i)**0.5) for i,w_ in enumerate(w)]
            vector_norm_distribution.append(w_c)
        elif found > 0:
            not_found += 1
        
        if not_found > min(1e4,len(db)): #if the interval between ''last found'' and ''not found'' more than 1e4 vectors, break
            break


    print(sample_vector_list)
    print("found %d vectors of length shorter than sqrt(4.0/3.0)*gh(L_%d). (expected %f)"%(found, rank, 1))
    print("min_sample_norm:%f"%(sqrt(min_vector_norm)))

    for i in range(len(vector_norm_distribution)):
        print("vector_norm_distribution[%d]: "%(i),vector_norm_distribution[i])
        l = sum(l_**2 for l_ in vector_norm_distribution[i])
        print("||w_c[%d]|| = %f"%(i,l**0.5))

        mean_value = np.average(vector_norm_distribution[i])
        std_value = np.std(vector_norm_distribution[i])
        print("mean_value of vector_norm_distribution[%d]: %f"%(i,mean_value))
        print("std_value of vector_norm_distribution[%d]: %f"%(i,std_value))

    T_FIND = time.time() - T1
    print("walltime: %f sec." % T_FIND)
    print("find finished.")
    print("\n\n")
'''
    
    '''
    tau = True
    pump_dim = 110
    while(tau):
        old_vector = g6k.M.get_r(0, 0)
        flast = pump(g6k, tracer, 0,pump_dim, 0,goal_r0=goal_r0, verbose = True)[0]
        if g6k.M.get_r(0, 0) < goal_r0  or pump_dim == 140:
            tau = False
            block_filename = "svpchallenge/pump-%02d-dim-%03d.txt" % (pump_dim,g6k.full_n)
            fn = open(block_filename, "w")
            fn.write(str(g6k.M.B))
            fn.close()
            break
        if g6k.M.get_r(0, 0) < old_vector:
            print("old vector's length = %.5f\n" %(old_vector**0.5))
            print("new vector's length = %.5f\n" %(g6k.M.get_r(0, 0)**0.5))
            block_filename = "svpchallenge/pump-%02d-dim-%03d.txt" % (pump_dim,g6k.full_n)
            fn = open(block_filename, "w")
            fn.write(str(g6k.M.B))
            fn.close()
            B = IntegerMatrix.from_file(block_filename)
            g6k = Siever(B, params, seed=seed)
        pump_dim = pump_dim + 10

    rank = 162
    flag = True
    while(flag):
        for i in range(1):
            flast = pump(g6k, tracer, 0,rank, round(rank/log(rank)),goal_r0=goal_r0, verbose = True)[0]
            print("---------current_initial_vector:%.5f--------"%(g6k.M.get_r(0, 0)**0.5))
            if g6k.M.get_r(0, 0)**0.5 < 165 or rank == 164:
                    flag = False
                    block_filename = "qarychallenge/pump-%02d-dim-%03d.txt" % (rank,g6k.full_n)
                    fn = open(block_filename, "w")
                    fn.write(str(g6k.M.B))
                    fn.close()
                    break
        rank = rank + 2

    rank = 139
    for i in range(3):
        flast = pump(g6k, tracer, 0,rank, round(rank/log(rank)),goal_r0=goal_r0, verbose = True)[0]
        print("---------current_initial_vector:%.5f--------"%(g6k.M.get_r(0, 0)**0.5))
        if g6k.M.get_r(0, 0)**0.5 < 137.56:
                block_filename = "qarychallenge/pump-%02d-dim-%03d.txt" % (rank,g6k.full_n)
                fn = open(block_filename, "w")
                fn.write(str(g6k.M.B))
                fn.close()
                break
    
    flast = pump(g6k, tracer, 0,183, 183-145,goal_r0=goal_r0, verbose = True)[0]
    block_filename = "qarychallenge/pump-%02d-dim-%03d.txt" % (183,330)
    fn = open(block_filename, "w")
    fn.write(str(g6k.M.B))
    fn.close()
    '''
    




    if verbose:
        logger.info("sol %d, %s" % (n, A[0]))

    norm = sum([x*x for x in A[0]])
    if verbose:
        logger.info("norm %.1f ,hf %.5f" % (norm**.5, (norm/gh_2)**.5))

    tracer.exit()


    if hasattr(tracer, "trace"):
        stat = tracer.trace
        stat.data["flast"] = flast
        return stat
    else:
        return None


def asvp():
    """
    Run a Workout until 1.05-approx-SVP on matrices with dimensions in ``range(lower_bound, upper_bound, step_size)``.
    """
    description = asvp.__doc__

    args, all_params = parse_args(description,
                                  load_matrix=None,
                                  pre_bkz=None,
                                  verbose=True,
                                  challenge_seed=0,
                                  workout__dim4free_dec=2,
                                  trace=True,
                                  gamma=1.40,          # 默认值
                                  preprocessed_file="svpchallenge/120-midmat-bkz49-pnjbkz-beta60-jump1-reduced-hf-1.428280-Tbkz-20.071112s-Tpnjbkz-751.049610s.txt")

    stats = run_all(asvp_kernel, all_params.values(),
                    lower_bound=args.lower_bound,
                    upper_bound=args.upper_bound,
                    step_size=args.step_size,
                    trials=args.trials,
                    workers=args.workers,
                    seed=args.seed)

    inverse_all_params = OrderedDict([(v, k) for (k, v) in all_params.items()])

    for (n, params) in stats:
        stat = stats[(n, params)]
        if stat[0] is None:
            logging.info("Trace disabled")
            continue

        if len(stat) > 0:
            cputime = sum([float(node["cputime"]) for node in stat])/len(stat)
            walltime = sum([float(node["walltime"]) for node in stat])/len(stat)
            flast = sum([float(node["flast"]) for node in stat])/len(stat)
            avr_db, max_db = db_stats(stat)
            fmt = "%48s :: m: %1d, n: %2d, cputime :%7.4fs, walltime :%7.4fs, flast : %2.2f, avr_max db: 2^%2.2f, max_max db: 2^%2.2f" # noqa
            print(fmt % (inverse_all_params[params], params.threads, n, cputime, walltime, flast, avr_db, max_db))
        else:
            logging.info("Trace disabled")

    if args.pickle:
        pickler.dump(stats, open("hkz-asvp-%d-%d-%d-%d.sobj" %
                                 (args.lower_bound, args.upper_bound, args.step_size, args.trials), "wb"))


if __name__ == '__main__':
    asvp()
