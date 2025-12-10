# -*- coding: utf-8 -*-
"""
BKZ Tours.
"""
import sys
from .pump import pump
#from .pump_cpu import pump
from .workout import workout
import six
import psutil  
import os
from math import log,pi,e,sqrt,exp
import time
from fpylll.util import gaussian_heuristic

from os import mkdir
# from ...pro_pnjBKZ_simulator.codes.util

try:
    basestring
except NameError:
    basestring = str


def get_current_slope(r, start_row=0, stop_row=-1):                         
    """
    A Python re-implementation of ``MatGSO.get_current_slope``.

        >>> from fpylll import IntegerMatrix, GSO, LLL, FPLLL
        >>> FPLLL.set_random_seed(1337)
        >>> A = IntegerMatrix.random(100, "qary", bits=30, k=50)
        >>> _ = LLL.reduction(A)
        >>> M = GSO.Mat(A); _ = M.update_gso()
        >>> from fpylll.tools.quality import get_current_slope
        >>> M.get_current_slope(0, 100)  # doctest: +ELLIPSIS
        -0.085500625...
        >>> get_current_slope(M.r(), 0, 100) # doctest: +ELLIPSIS
        -0.085500625...

    """
    x = [log(r[i]) for i in range(start_row, stop_row)]                     
    n = stop_row - start_row                                                
    i_mean = (n - 1) * 0.5 + start_row                                      
    x_mean = sum(x)/n                                                       
    v1, v2 = 0.0, 0.0
    for i in range(stop_row - start_row):
        v1 += (i - i_mean) * (x[i] - x_mean)
        v2 += (i - i_mean) * (i - i_mean)
    return v1 / v2                                                          

def dim4free_wrapper(dim4free_fun, blocksize):
    """
    Deals with correct dim4free choices for edge cases when non default
    function is chosen.

    :param dim4free_fun: the function for choosing the amount of dim4free
    :param blocksize: the BKZ blocksize

    """
    if blocksize < 40:
        return 0
    dim4free = dim4free_fun(blocksize)
    return int(min((blocksize - 40)/2, dim4free))


def default_dim4free_fun(blocksize):
    """
    Return expected number of dimensions for free, from exact-SVP experiments.

    :param blocksize: the BKZ blocksize

    """
    return int(11.5 + 0.075*blocksize)

def theo_dim4free_fun(blocksize):
     """
     Theoretical Dimension-for-free function in [Duc18]
     """

     return int(blocksize*log(4/3.)/log(blocksize/(2.*pi*e))) 


def theo_dim4free_fun1(blocksize):
     """
     Theoretical Dimension-for-free function in [Duc18]
     """

     return int(blocksize*log(4/3.)/log(blocksize/(2.*pi))) 


def theo_dim4free_fun2(blocksize,rhf):
     """
     wlz asiaccs23 pessimic d4f
     """

     for f in range(1,int(blocksize/2)):
            right = log(sqrt((4*(blocksize-f))/(3*blocksize)))/log(rhf)
            if f >=right:
                return f


def theo_dim4free_fun3(blocksize,rhf):
     """
      wlz asiaccs23 optimistic d4f
     """

     for f in range(1,int(blocksize/2)):
            right = log(sqrt(4.0/3.0))/log(rhf)
            if f >=right:
                return f


#wnjbkz&headwnjbkz
def workout_bkz_tour(g6k, tracer,blocksize, dim4free_fun=theo_dim4free_fun2,rank=0,
                   extra_dim4free=0, jump = 2, tours = 1, total_tours = 1,workout_params=None, pump_params=None, goal_r0=0.):
    """
    Run a workout BKZ-tour: call ``workout`` as an SVP oracle consecutively on
    each block.

    :param g6k: The g6k object to work with
    :param tracer: A tracer for g6k
    :param blocksize: dimension of the blocks
    :param dim4free_fun: number of dimension for free as a function of beta (function, or string e.g. `lambda x: 11.5+0.075*x`)
    :param extra_dim4free: increase the number of dims 4 free (blocksize is increased, but not sieve dimension)
    :param workout_params: parameters to pass to the workout
    :param pump_params: parameters to pass to the pump

    """
    if workout_params is None:
        workout_params = {}


    if "dim4free_min" in workout_params:
        raise ValueError("In naive_bkz, you should choose dim4free via dim4free_fun.")

    d = g6k.full_n

    if isinstance(dim4free_fun, basestring):
        dim4free_fun = eval(dim4free_fun)
    

    r = [g6k.M.get_r(j, j) for j in range(d)]
    log_vol = sum([log(x) for x in r])
    log_vol_n = 1.0 / d * log_vol
    vol_n = sqrt(exp(log_vol_n))
    hf = g6k.M.get_r(0,0)**0.5/vol_n
    rhf = pow(hf, 1.0/d)

    dim4free_min = theo_dim4free_fun2(blocksize,rhf) + extra_dim4free
    dim4free_max = theo_dim4free_fun3(blocksize,rhf) + extra_dim4free
    print(extra_dim4free)
    print(dim4free_min)
    print(dim4free_max)

    inner_tours = jump

    blocksize += extra_dim4free
    jump_kappa = jump

    if rank ==0:
        if tours == 0:
            kappa_list =list(range(0,int((d-50-dim4free_max)/jump_kappa)*jump_kappa+1,jump_kappa)) # kappa must end on at least the d - 50 -f_max (f_max is in the workout)
        else:
            kappa_list =list(range(0,int((d-blocksize+dim4free_max)/jump_kappa)*jump_kappa+1,jump_kappa)) # kappa must end on at least the d - 50 -f_max (f_max is in the workout)
    else:
        if tours == 0:
            kappa_list =list(range(0,int(min(rank,d-50-dim4free_max)/jump_kappa)*jump_kappa+1,jump_kappa)) # kappa must end on at least the d - 50 -f_max (f_max is in the workout)
        else:
            kappa_list =list(range(0,int(min(rank,d-blocksize+dim4free_max)/jump_kappa)*jump_kappa+1,jump_kappa)) # kappa must end on at least the d - 50 -f_max (f_max is in the workout)

    print(kappa_list)
    flast = 100
    for kappa in kappa_list:
        beta = min(blocksize, d- kappa)
        if beta >=blocksize:
            f_min = min(dim4free_min,d-kappa)
            f_max= min(dim4free_max,d-kappa)
            print(f_min)
            print(f_max)
        else:
            f_min = theo_dim4free_fun2(beta,rhf) + extra_dim4free
            f_max = theo_dim4free_fun3(beta,rhf) + extra_dim4free
            print(f_min)
            print(f_max)
            inner_tours = jump
        for t in range(inner_tours):
            print("\r %dth inner_tours k:%d, b:%d, f_max:%d , RAM cost: %.4f GB" % (t,kappa, beta, f_max, psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024), end=' ')
            sys.stdout.flush()
            flast = min(flast,workout(g6k, tracer, rank,kappa, beta, f_min,f_max,goal_r0=goal_r0,pump_params=pump_params, dim4free_dec=jump,current_beta = blocksize,extra_dim4free=extra_dim4free))
            g6k.lll(0, d)
    return flast
       
        


def pump_n_jump_bkz_tour(g6k, tracer, blocksize,  jump=1,
                         dim4free_fun=theo_dim4free_fun, extra_dim4free=0,
                         pump_params=None, goal_r0=0., verbose=False):
    """
    Run a PumpNjump BKZ-tour: call Pump consecutively on every (jth) block.

    :param g6k: The g6k object to work with
    :param tracer: A tracer for g6k
    :param blocksize: dimension of the blocks
    :param jump: only call the pump every j blocks
    :param dim4free_fun: number of dimension for free as a function of beta (function, or string
        e.g. `lambda x: 11.5+0.075*x`)
    :param extra_dim4free: increase the number of dims 4 free (blocksize is increased, but not sieve
        dimension)
    :param pump_params: parameters to pass to the pump
    """
    if pump_params is None:
        pump_params = {"down_sieve": False}

    if "dim4free" in pump_params:
        raise ValueError("In pump_n_jump_bkz, you should choose dim4free via dim4free_fun.")

    d = g6k.full_n
    g6k.shrink_db(0)
    g6k.lll(0,d)
    g6k.update_gso(0,d)

    if isinstance(dim4free_fun, six.string_types):
        dim4free_fun = eval(dim4free_fun)
    
    bkz_betas = [blocksize,jump]
    r = [g6k.M.get_r(j, j) for j in range(d)]
    log_vol = sum([log(x) for x in r])
    log_vol_n = 1.0 / d * log_vol
    vol_n = sqrt(exp(log_vol_n))
    hf = g6k.M.get_r(0,0)**0.5/vol_n
    rhf = pow(hf, 1.0/d)

    
    #original_pnj_bkz
    dim4free = dim4free_wrapper(dim4free_fun, blocksize) + extra_dim4free
    blocksize += extra_dim4free
    indices  = [(0, blocksize - dim4free + i, i) for i in range(0, dim4free, jump)]
    indices += [(i, blocksize, dim4free) for i in range(0, d - blocksize, jump)]
    indices += [(d - blocksize + i, blocksize - i, dim4free - i) for i in range(0, dim4free, jump)]
    print(indices)


    pump_params["down_stop"] = dim4free+3


    max_RAM_cost = 0
    

    for (kappa, beta, f) in indices:
        if verbose:
            print("\r k:%d, b:%d, f:%d , RAM cost: %.4f GB" % (kappa, beta, f, psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024), end=' ')
            sys.stdout.flush()
            
            RAM_cost = psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024
            if max_RAM_cost < RAM_cost:
                max_RAM_cost = RAM_cost
        
        flast = pump(g6k, tracer, kappa, beta, f, **pump_params)[0]

        sub_dir = "qarychallenge/dim-%d/pnjbkz-%d-jump-%d" % (d, blocksize - extra_dim4free, jump)
        os.makedirs(sub_dir, exist_ok=True)

        if kappa >= 0:
            block_filename = "qarychallenge/dim-%d/pnjbkz-%d-jump-%d/pnjbkz-%d-jump-%d-kappa-%d-beta-%d-f-%d" % (d,blocksize - extra_dim4free, jump,blocksize - extra_dim4free, jump,kappa,beta,f)
            fn = open(block_filename, "w")
            fn.write(str(g6k.M.B))
            fn.close()
        
        g6k.lll(0, d)
        if g6k.M.get_r(0, 0) <= goal_r0:
            return max_RAM_cost

    if verbose:
        
        print("\r k:%d, b:%d, f:%d  , RAM cost: %.4f GB " % (d-(blocksize-dim4free), blocksize-dim4free, 0,psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024), end=' ')
        sys.stdout.flush()

        RAM_cost = psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024
        if max_RAM_cost < RAM_cost:
            max_RAM_cost = RAM_cost

    

    if verbose:
        print("\r k:%d, b:%d, f:%d " % (d-(blocksize-dim4free), blocksize-dim4free, 0))
        
        sys.stdout.flush()

    pump_params["down_stop"] = blocksize - dim4free
    
    T_0 = time.time()
    pump(g6k, tracer, d-(blocksize-dim4free), blocksize-dim4free, 0, **pump_params)
       
 
    if verbose:
        print('')
        sys.stdout.flush()
    return flast,max_RAM_cost

