from math import *
from model_BKZ import *
from proba_util import gaussian_center_weight

log_infinity = 9999

STEPS_b = 1
STEPS_m = 5

from math import log, exp, lgamma, pi

from decimal import Decimal, getcontext

def calculate_target_log(L_b, L_s, precision=50):
    """
    Calculate target log2(1/(p_sieve*(1-p_bkz)+p_bkz)) based on known log2(1/p_bkz) and log2(1/p_sieve)
    
    Parameters:
        L_b (float/Decimal): log2(1.0 / p_bkz), known value
        L_s (float/Decimal): log2(1.0 / p_sieve), known value
        precision (int): Calculation precision (number of significant digits for Decimal module, default 50, adjustable as needed)
    
    Returns:
        Decimal: Target logarithmic result log2(1/(p_sieve*(1-p_bkz)+p_bkz))
    """
    # 1. Configure Decimal precision (avoid precision loss when calculating small probabilities)
    getcontext().prec = precision
    # Extend maximum exponent range (prevent overflow when 2^L_b is too large)
    getcontext().Emax = 999999

    # 2. Convert input L_b, L_s to Decimal type (ensure consistent operation types)
    # If input is float, first convert to string then to Decimal to avoid floating point precision loss
    if isinstance(L_b, float):
        L_b = Decimal(str(L_b))
    if isinstance(L_s, float):
        L_s = Decimal(str(L_s))

    # 3. Calculate original probabilities p_bkz and p_sieve
    # p_bkz = 1 / (2^L_b), p_sieve = 1 / (2^L_s)
    p_bkz = Decimal(1) / (Decimal(2) ** L_b)
    p_sieve = Decimal(1) / (Decimal(2) ** L_s)

    # 4. Calculate inner denominator D = p_sieve*(1 - p_bkz) + p_bkz of target expression
    term1 = p_sieve * (Decimal(1) - p_bkz)  # p_sieve*(1-p_bkz)
    D = term1 + p_bkz                       # Final denominator D

    # 5. Calculate target logarithm: log2(1/D) = -log2(D)
    # Use Decimal.log10() to indirectly calculate log2 (since Decimal has no direct log2 method)
    # Formula: log2(x) = log10(x) / log10(2)
    log2_D = D.log10() / Decimal(2).log10()
    L_target = -log2_D  # Target result

    return L_target


def ball_log_vol(n):
    """
    Return volume of `n`-dimensional unit ball

    :param n: dimension

    """
    return (n/2.) * log(pi) - lgamma(n/2. + 1)

def gaussian_heuristic(q,h,w):
    """
    Return squared norm of shortest vector as predicted by the Gaussian heuristic.

    :param r: vector of squared Gram-Schmidt norms

    """
    log_vol = h*log(q)
    log_gh =  1./w * (log_vol - ball_log_vol(w))
    return exp(log_gh)



def theo_dim4free_fun_optimistic(blocksize):
     """
     Theoretical Dimension-for-free function in [Duc18]
     """
     return float(blocksize*log(4/3.)/log(blocksize/2./pi/e)) 

def compare_exponents(a, b, c, d):
    # Use maximum exponent for scaling instead of minimum exponent
    max1 = max(a, b)
    max2 = max(c, d)
    
    # Calculate scaled terms to avoid direct calculation of large exponents
    term1 = pow(2, a - max1) + pow(2, b - max1)
    term2 = pow(2, c - max2) + pow(2, d - max2)
    
    # Calculate sum in logarithmic form
    sum1_log = max1 + log2(term1)
    sum2_log = max2 + log2(term2)
    
    return sum1_log > sum2_log


class MSISParameterSet:
    def __init__(self, n, w, h, B, q, norm=""):
        self.n = n      # Ring Dimension
        self.w = w      # MSIS dimension
        self.h = h      # Number of equations
        self.B = B      # norm bound
        self.q = q      # Modulus
        self.norm = norm


def SIS_l2_cost(q, w, h, B, b, cost_svp=svp_classical, verbose=False):
    """ Return the cost of finding a vector shorter than B with BKZ-b if it works.
    The equation is Ax = 0 mod q, where A has h rows, and w collumns (h equations in dim w).
    """
    if B>=q:
        if verbose:
            print ("Do not know how to handle B>=q in l2 norm. Concluding 0 bits of security.")
        return 0
    l = BKZ_first_length(q, h, w-h, b)
    if l > B:
        return log_infinity 
    if verbose:
        print ("Attack uses block-size %d and %d equations"%(b, h))
        print ("shortest vector used has length l=%.2f, q=%d, `l<q'= %d"%(l, q, l<q))
    return cost_svp(b)

def SIS_linf_cost_bai(q, w, h, B, b, cost_svp=svp_classical, verbose=False):
    """ Return the cost of finding a vector shorter than B in infinity norm, using BKZ-b, if it works.
    The equation is Ax = 0 mod q, where A has h rows, and w columns (h equations in dim w).
    """
    
    (i, j, L) = construct_BKZ_shape_randomized(q, h, w-h, b)
    #(i, j, L) = construct_BKZ_shape(h * log(q), 1, w-1, b)

    #print("(i,j):(%d,%d)"%(i,j))

    l = exp(L[i])

    d = j - i + 1
    #sigma = l / sqrt(j - i + 1)
    #sigma = (sqrt(4.0/3.0)*l)/ sqrt(j - i + 1)  
    sigma = l/ sqrt(j - i + 1)
    p_middle = gaussian_center_weight(sigma, B)
    p_head = 2.*B / q

    log2_eps = d * log(p_middle, 2) + i * log(p_head, 2)
    log2_R = max(0, - log2_eps - nvec_sieve(b))

    if verbose:
        print ("Attack uses block-size %d and %d dimensions, with %d q-vectors"%(b, w, i))
        print ("log2(epsilon) = %.2f, log2 nvector per run %.2f"%(log2_eps, nvec_sieve(b)))
        print ("shortest vector used has length l=%.2f, q=%d, `l<q'= %d"%(l, q, l<q))
    return cost_svp(b) + log2_R, b


def SIS_linf_cost_sub(i,j,l,B,q,b,C,cost_svp=svp_classical):
    """ Return the cost of finding a vector shorter than B in infinity norm, using BKZ-b, if it works.
    The equation is Ax = 0 mod q, where A has h rows, and w columns (h equations in dim w).
    """

    d = j - i + 1

    sigma = l/ sqrt(j - i + 1)
    p_middle = gaussian_center_weight(sigma, B)
    p_head = 2.*B / q

    log2_eps = d * log(p_middle, 2) + i * log(p_head, 2)
    log2_P = max(0, -log2_eps - nvec_sieve(b))
    log2_C = log(C,2)

    return cost_svp(b) + log2_C, log2_P


def SIS_linf_cost_our(q, w, h, B, b, cost_svp=svp_classical, verbose=False):
    """ Return the cost of finding a vector shorter than B in infinity norm, using BKZ-b, if it works.
    The equation is Ax = 0 mod q, where A has h rows, and w columns (h equations in dim w).
    """
    
    (i, j, L) = construct_BKZ_shape_randomized(q, h, w-h, b)

    l = exp(L[i])

    rank = 9999
    gh = gaussian_heuristic(q,h,w)

    alpha = l/gh
    for i_prime in range(50,w,1):
        right = pow(8*log(10),1.0/b)*sqrt(i_prime/w) * pow(delta_BKZ(b),w - i_prime + 1) 
        if (alpha >= right):
            rank = i_prime
            break
        
    bkz_cost,p_bkz = SIS_linf_cost_sub(i,j,l,B,q,b,1)
    sieve_cost,p_sieve = SIS_linf_cost_sub(i,j,l,B,q,rank,1)

    # To enable D4f as used in the paper, uncomment these lines:
    #bkz_cost, p_bkz = SIS_linf_cost_sub(i,j,l,B,q,round(b - theo_dim4free_fun_optimistic(b)),1)
    #sieve_cost, p_sieve = SIS_linf_cost_sub(i,j,l,B,q,round(rank - theo_dim4free_fun_optimistic(rank)),1)

    bkz_cost_dec = Decimal(str(bkz_cost)) if isinstance(bkz_cost, float) else bkz_cost
    sieve_cost_dec = Decimal(str(sieve_cost)) if isinstance(sieve_cost, float) else sieve_cost
    
    # Calculate max1 (now Decimal type)
    max1 = max(bkz_cost_dec, sieve_cost_dec)

   # Now both operands of subtraction operation are Decimal type
    term1 = pow(2, bkz_cost_dec - max1) + pow(2, sieve_cost_dec - max1)
    
    # Calculate log2_p_total (unchanged)
    log2_p_total = calculate_target_log(p_bkz, p_sieve, precision=50)
    
    # Calculate log2(term1) and convert to Decimal type
    log2_term1 = Decimal(log(term1)) / Decimal(log(2))
    
    # All variables are Decimal type, can be added normally
    Total_cost = log2_p_total + max1 + log2_term1

    flag = True
    for b_prime in range(b+1,w):
        (i_prime, j_prime, L_prime) = construct_BKZ_shape_randomized(q, h, w-h, b_prime)

        l_prime = exp(L_prime[i])

        rank_prime = 9999

        alpha_prime = l_prime/gh
        for i in range(50,w,1):
            right = pow(8*log(10),1.0/b_prime)*sqrt(i/w) * pow(delta_BKZ(b_prime),w - i + 1) 
            if (alpha_prime >= right):
                rank_prime = i
                break
    
        bkz_cost_prime,p_bkz_prime = SIS_linf_cost_sub(i_prime,j_prime,l_prime,B,q,b_prime,1)
        sieve_cost_prime,p_sieve_prime = SIS_linf_cost_sub(i_prime,j_prime,l_prime,B,q,rank_prime,1)
        
        # To enable D4f as used in the paper, uncomment these lines:
        #bkz_cost_prime, p_bkz_prime = SIS_linf_cost_sub(i_prime,j_prime,l_prime,B,q,round(b_prime - theo_dim4free_fun_optimistic(b_prime)),1)
        #sieve_cost_prime, p_sieve_prime = SIS_linf_cost_sub(i_prime,j_prime,l_prime,B,q,round(rank_prime - theo_dim4free_fun_optimistic(rank_prime)),1)
        
        
        bkz_cost_dec_prime = Decimal(str(bkz_cost_prime)) if isinstance(bkz_cost_prime, float) else bkz_cost_prime
        sieve_cost_dec_prime = Decimal(str(sieve_cost_prime)) if isinstance(sieve_cost_prime, float) else sieve_cost_prime
        
        max2 = max(bkz_cost_dec_prime, sieve_cost_dec_prime)

    # Now both operands of subtraction operation are Decimal type
        term2 = pow(2, bkz_cost_dec_prime - max2) + pow(2, sieve_cost_dec_prime - max2)
        
        # Calculate log2_p_total (unchanged)
        log2_p_total_prime = calculate_target_log(p_bkz_prime, p_sieve_prime, precision=50)
        
        # Calculate log2(term1) and convert to Decimal type
        log2_term2 = Decimal(log(term2)) / Decimal(log(2))
        
        # All variables are Decimal type, can be added normally
        Total_cost_prime = log2_p_total_prime + max2 + log2_term2
        #if compare_exponents(bkz_cost, sieve_cost, bkz_cost_prime, sieve_cost_prime):
        if Total_cost > Total_cost_prime:
            flag = False
            break
    
    if flag == True:
        return Total_cost,max(b,rank),flag
        #return Total_cost,max(round(b - theo_dim4free_fun_optimistic(b)),round(rank - theo_dim4free_fun_optimistic(rank))),flag
    else:
        return log_infinity,50,flag
        
    

    if verbose:
        print ("Attack uses block-size %d and %d dimensions, with %d q-vectors"%(b, w, i))
        print ("log2(epsilon) = %.2f, log2 nvector per run %.2f"%(log2_eps, nvec_sieve(b)))
        print ("shortest vector used has length l=%.2f, q=%d, `l<q'= %d"%(l, q, l<q))
    return cost_svp(b) + log2_R, b

def SIS_optimize_attack_our(q, max_w, h, B, cost_attack=SIS_linf_cost_our, cost_svp=svp_classical, verbose=False):
    """ Find optimal parameters for a given attack
    """
    best_cost = log_infinity

    for b in range(50+30, max_w, STEPS_b):
        for w in [max_w]:    # No need to exhaust w here as the attack will auto-adjust anyway  range(max(h+1, b+1), max_w, STEPS_m):
            best_cost, best_b, flag = cost_attack(q, w, h, B, b, cost_svp)
            if flag:
                return (w, best_b, best_cost)


def SIS_optimize_attack(q, max_w, h, B, cost_attack=SIS_linf_cost_bai, cost_svp=svp_classical, verbose=False):
    """ Find optimal parameters for a given attack
    """
    best_cost = log_infinity

    for b in range(50, max_w, STEPS_b):
        if cost_svp(b) > best_cost:
            break
        for w in [max_w]:    # No need to exhaust w here as the attack will auto-adjust anyway  range(max(h+1, b+1), max_w, STEPS_m):
            cost, rank = cost_attack(q, w, h, B, b, cost_svp)
            if cost<=best_cost:
                best_cost = cost
                best_w = w
                best_b = rank

    if verbose:
        cost_attack(q, best_w, h, B, best_b, cost_svp=cost_svp, verbose=verbose)

    return (best_w, best_b, best_cost)



def check_eq(m_pc, m_pq, m_pp):
    if (m_pc != m_pq):
        print("m and b not equals among the three models")
    if (m_pq != m_pp):
        print("m and b not equals among the three models")


def MSIS_summarize_attacks(ps):
    """ Create a report on the best primal and dual BKZ attacks on an l_oo - MSIS instance
    """
    q = ps.q
    h = ps.n * ps.h
    max_w = ps.n * ps.w
    B = ps.B

    if ps.norm == "linf":
        attack = SIS_linf_cost_our
        #attack = SIS_linf_cost_bai
    elif ps.norm == "l2":
        attack = SIS_l2_cost
    else:
        raise ValueError("Unknown norm: "+ps.norm)

    (m_pc, b_pc, c_pc) = SIS_optimize_attack_our(q, max_w, h, B, cost_attack=attack, cost_svp=svp_classical, verbose=False)
    (m_pq, b_pq, c_pq) = SIS_optimize_attack_our(q, max_w, h, B, cost_attack=attack, cost_svp=svp_quantum, verbose=False)
    (m_pp, b_pp, c_pp) = SIS_optimize_attack_our(q, max_w, h, B, cost_attack=attack, cost_svp=svp_plausible, verbose=False)

    print("b_pc:%d"%b_pc)
    print("b_pc:%d"%b_pq)
    print("b_pp:%d"%b_pp)

    check_eq(m_pc, m_pq, m_pp)
    check_eq(b_pc, b_pq, b_pp)

    #print("SIS & %d & %d & %d & %d & %d"%(m_pq, b_pq, int(floor(c_pc)), int(floor(c_pq)), int(floor(c_pp))))
    print("SIS & %d & %d & %.2f & %.2f & %.2f"%(m_pq, b_pq, c_pc, c_pq, c_pp))

    #return (b_pq, int(floor(c_pc)), int(floor(c_pq)), int(floor(c_pp)))
    return (b_pq, c_pc, c_pq, c_pp)