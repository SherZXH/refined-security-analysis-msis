# -*- coding: utf-8 -*-
"""

"""
import sys
from .pump import pump
from fpylll.util import gaussian_heuristic
import time
from math import log,pi ,e
import os

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


def workout(g6k, tracer, kappa, blocksize, dim4free_min=0,              # Main parameters
            dim4free_max=0,dim4free_dec=2,current_beta = 50, extra_dim4free=0,start_n=50, goal_r0=0.,                     # Loop control
            verbose= False, save_prefix=None, prefer_left_insert=1.04,pump_params=None           # Misc
            ):
    """
    :param g6k: The g6k object to work with
    :param tracer: A tracer for g6k
    :param kappa: beginning of the block
    :param blocksize: dimension of the block
    :param dim4free_min: Minimal number of dimension for free ``dimension for free'' [Ducas,
        Eurcrypt 2018] (may stop before reaching that if goal_r0)
    :param dim4free_dec: By how much do we decreaseee dim4free at each iteration
    :param start_n: Dimension of the first pump
    :param goal_r0: an extra hook to always insert at position kappa if this goal length can be met
        by a lift.  Quit when this is reached.
    :param verbose: Print workout steps (with timing and quality) information on the standard
        output.  Enforce verbosity of pump as well.
    :param save_prefix: If not None, save intermediate basis at a file-name with this prefix.
        Allows to resume computation.
    :param pump_params: Para fs = list(range(dim4free_min, f_start, dim4free_dec)[::-1])    meters to forward to the pump.

    """
    if pump_params is None:
        pump_params = {}
    '''
    f_start = max(blocksize - start_n, 0, dim4free_min)
    '''
    
    f_start = 0
    print(f_start)
    if kappa == 0:
        fs = list(range(f_start, dim4free_max+1,dim4free_dec))    
    else:
        fs = list(range(dim4free_max, dim4free_min,-dim4free_dec))
        if fs[-1] !=  dim4free_min:
            fs.append(dim4free_min)
    
    print(fs)
    '''
    fs = list(range(dim4free_min, f_start+1, -round(dim4free_min/6)))           # could not arrive the target saturation-ratio 0.5 and the practical saturation-ratio is 0.4x from 192
    print(fs)
    '''
    
    """
    if goal_r0:
        fs += 3*[dim4free_min]                                   
    """
    gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(kappa, kappa+blocksize)])
    runtimestart = time.time()                                    
                                                                            
    if "verbose" not in pump_params:
        pump_params["verbose"] = verbose                            
    pump_params["prefer_left_insert"] = prefer_left_insert

    flast = 100
    with tracer.context(("workout", "kappa:%d beta:%d f:%d" % (kappa, blocksize, dim4free_max))):
        for f in fs:
            timestart = time.time()                                 # 每轮计时

            sys.stdout.flush()                                      # 把当前处理过程打印到屏幕上
            if kappa ==0:
                print("\r workout-   k:%d, blocksize:%d, fs:%d " % (kappa, blocksize-dim4free_max+f,f))
            else:
                print("\r workout-   k:%d, blocksize:%d, fs:%d " % (kappa, blocksize,f))

            if kappa ==0:
                flast = min(flast, pump(g6k, tracer, kappa, blocksize-dim4free_max+f, f, goal_r0=goal_r0, down_sieve =  True,verbose = False)[0])      # 调用PUMP算法 
            else:
                flast = min(flast, pump(g6k, tracer, kappa, blocksize, f, goal_r0=goal_r0, down_sieve = True, verbose = False)[0])      # 调用PUMP算法 
            '''
            if not os.path.isdir("svpchallenge/dim-%d/workoutbkz-%d-jump-%d" % (g6k.full_n,blocksize - extra_dim4free,dim4free_dec)):
                os.mkdir("svpchallenge/dim-%d/workoutbkz-%d-jump-%d" % (g6k.full_n,blocksize - extra_dim4free,dim4free_dec))
            '''
            
            if not os.path.isdir("qarychallenge/dim-%d/workoutbkz-%d-jump-%d" % (g6k.full_n,blocksize - extra_dim4free,dim4free_dec)):
                os.mkdir("qarychallenge/dim-%d/workoutbkz-%d-jump-%d" % (g6k.full_n,blocksize - extra_dim4free,dim4free_dec))

            if kappa == 0:
                #block_filename = "svpchallenge/dim-%d/workoutbkz-%d-jump-%d/workoutbkz-%d-jump-%d-kappa-%d-beta-%d-f-%d" % (g6k.full_n,current_beta - extra_dim4free,dim4free_dec,blocksize - extra_dim4free,dim4free_dec,kappa,blocksize-dim4free_max+f,f)
                block_filename = "qarychallenge/dim-%d/workoutbkz-%d-jump-%d/workoutbkz-%d-jump-%d-kappa-%d-beta-%d-f-%d" % (g6k.full_n,current_beta - extra_dim4free,dim4free_dec,blocksize - extra_dim4free,dim4free_dec,kappa,blocksize-dim4free_max+f,f)
                fn = open(block_filename, "w")
                fn.write(str(g6k.M.B))
                fn.close()
            else:
                #block_filename = "svpchallenge/dim-%d/workoutbkz-%d-jump-%d/workoutbkz-%d-jump-%d-kappa-%d-beta-%d-f-%d" % (g6k.full_n,current_beta - extra_dim4free,dim4free_dec,blocksize - extra_dim4free,dim4free_dec,kappa,blocksize,f)
                block_filename = "qarychallenge/dim-%d/workoutbkz-%d-jump-%d/workoutbkz-%d-jump-%d-kappa-%d-beta-%d-f-%d" % (g6k.full_n,current_beta - extra_dim4free,dim4free_dec,blocksize - extra_dim4free,dim4free_dec,kappa,blocksize,f)
                fn = open(block_filename, "w")
                fn.write(str(g6k.M.B))
                fn.close()

            """
            native_bkz_filename = "qarychallenge/dim-330/gpu/nativebkz/blocksize-%d/nativebkz-kappa-%d-beta-%d-f-%d" % (kappa, blocksize, dim4free_min) #store after native-bkz processing lattice basis matrix step by step
            fn = open(native_bkz_filename, "w")
            fn.write(str(g6k.M.B))
            fn.close()
             """
            if verbose:
                gh2 = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(kappa+f, kappa+blocksize)])           
                quality = (gh * (blocksize - f)) / (gh2 * blocksize)                                             
                r = (gh) / (gh2)
                print("T:%10.5fs, TT:%10.5fs, q:%10.5f r:%10.5f r0/gh:%10.5f" %
                      (time.time() - timestart,
                       time.time() - runtimestart, quality,r, g6k.M.get_r(kappa, kappa) / gh))                    
                #if f == dim4free_min:
                #   fn = open("experiments_dh_hash/%d_%d_%f_%d_%d_%d.txt" % (dim4free_min,blocksize,prefer_left_insert,g6k.params.dh_dim,g6k.params.dh_min), "w")                      
                #   fn.write(str(r))
                #   fn.close()
                #print >> sys.stderr, quality
                #print >> sys.stderr, g6k.M.B

            if g6k.M.get_r(0, 0) < goal_r0:                                                               
                if save_prefix is not None:
                    fn = open("%s_%d_%d.sol" % (save_prefix.rstrip(), g6k.M.d - f, g6k.M.d), "w")                 
                    fn.write(str(g6k.M.B))
                    fn.close()
                break

            if save_prefix is not None:
                fn = open("%s_%d_%d.mat" % (save_prefix.rstrip(), g6k.M.d - f, g6k.M.d), "w")                      
                fn.write(str(g6k.M.B))
                fn.close()

    return flast
