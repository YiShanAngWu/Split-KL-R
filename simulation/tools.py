# 24 Jan, 2022
# Some helper functions

import numpy as np
from scipy import optimize
from math import log, sqrt, ceil, exp
from scipy.stats import binom

"""
KL Inequalities
"""

def KL(Q, P):
    """
    Compute Kullback-Leibler (KL) divergence between distributions Q and P.
    """
    return sum([ q*log(q/p) if q > 0. else 0. for q,p in zip(Q,P) ])
    
def KL_binomial(q, p):
    """
    Compute the KL-divergence between two Bernoulli distributions of probability
    of success q and p. That is, Q=(q,1-q), P=(p,1-p).
    """
    return KL([q, 1.-q], [p, 1.-p])
    
def solve_kl_sup(q, right_hand_side):
    """
    find x such that:
        kl( q || x ) = right_hand_side
        x > q
    """
    f = lambda x: KL_binomial(q, x) - right_hand_side

    if f(1.0-1e-9) <= 0.0:
        return 1.0-1e-9
    else:
        return optimize.brentq(f, q, 1.0-1e-9)

def solve_kl_inf(q, right_hand_side):
    """
    find x such that:
        kl( q || x ) = right_hand_side
        x < q
    """
    f = lambda x: KL_binomial(q, x) - right_hand_side

    if f(1e-9) <= 0.0:
        return 1e-9
    else:
        return optimize.brentq(f, 1e-9, q)

"""
Inverse of Chernoff-Hoeffding's Bound
"""
def solve_CH_left(hatp, right_hand_side):
    """
    find p such that:
        eps = p - hatp >= 0
        kl( p-eps || p) = right_hand_side
        eps < p
    """
    f = lambda x: KL_binomial(hatp, x) - right_hand_side

    if f(1-1e-9) >= 0.0:
        return 1-1e-9
    else:
        return optimize.brentq(f, 1e-9, 1)

def solve_CH_right(p, right_hand_side):
    """
    find eps such that:
        kl( p+eps || p) = right_hand_side
        eps < 1-p
    """
    f = lambda x: KL_binomial(p+x, p) - right_hand_side

    if f(1-p-1e-9) >= 0.0:
        return 1-p-1e-9
    else:
        return optimize.brentq(f, 1e-9, 1-p)

"""
Binomial Tail
"""
    
def solve_BinL(n, k, delta):
    """
    find p such that:
        binom.cdf(k, n, p) > delta
        k/n < p <= 1
    """
    f = lambda x: binom.cdf(k, n, x) - delta
    
    if f(1.0) >= 0.0:
        return 1.0
    else:
        return optimize.brentq(f, k/n, 1.0)
        
def solve_BinR(n, k, delta):
    """
    find p such that:
        1 - binom.cdf(k-1, n, p) > delta
        0 <= p < k/n
    """
    f = lambda x: 1-binom.cdf(k-1, n, x) - delta
    
    if k ==0:
        return 0.0
    elif f(0.0) >= 0.0:
        return 0.0
    else:
        return optimize.brentq(f, 0., k/n)

"""
Bennett's Inequalities
"""

def OBennett(q, Var, n, delta):
    """
    Bennett's inequality with oracle variance Var
    """
    return min(1-1e-9, q + sqrt(2*Var*log(1/delta)/n) + 2*log(1/delta)/(3*n))


def Emp_Bernstein(q, hVar, a, b, n, delta):
    """
    Empirical Bennett's inequality (with empirical variance hVar)
    generalize to Z\in[a,b]
    """
    return min(b, q + sqrt(2*hVar*log(2/delta)/n) + 7*(b-a)*log(2/delta)/(3*n))

def Unexp_Bernstein(q, x2, b, n, delta):
    """
    Unexpected Bernstein's inequality
    """
    gamGridSize = ceil(log(0.5*sqrt(n/log(1/delta)))/log(2))
    bounds = np.zeros(gamGridSize)
    for i in range(gamGridSize):
        gam = 1/(b*2**(i+1)) 
        bounds[i] = q + (-log(1-gam*b)/(gam*b*b) - 1/b)*x2 + (log(gamGridSize/delta))/(gam*n)
    return min(b, min(bounds))

"""
Empirical Mean and Variance
"""

def Emp_MaV(vec):
    """
    Compute the empirical mean and the empirical variance
    of a sequence of random variables vec
    """
    mean = np.mean(vec)
    var = np.sum((vec-mean)**2)/len(vec-1)
    x2 = np.mean(vec**2)
    return mean, var, x2
    
"""
kl shifted bound
"""
def splitkl(mu, a, b, barSzp, barSzm, right_hand_side):
    """
    Compute the kl shifted bound.
    mu = shift
    Z=(z1,z2,z3), the support of the ternary random variable
    barSzp, barSzm = \hat{Z}_n^+, \hat{Z}_n^-
    """
    barSzpT = 0 if mu==b else barSzp/(b-mu)
    barSzmT = 0 if mu==a else barSzm/(mu-a)
    EzpSUP = solve_kl_sup(barSzpT, right_hand_side)
    EzmINF = solve_kl_inf(barSzmT, right_hand_side)
    return mu + (b-mu)*EzpSUP - (mu-a)*EzmINF

"""
Catoni's bound
"""
def Catoni(q, n, delta):
    CGridSize = ceil(log(0.5*sqrt(n/log(1/delta)))/log(2))
    eps = log(CGridSize/delta)/n
    bounds = np.zeros(CGridSize)
    for i in range(CGridSize):
        C = 1/2**i
        bounds[i] = (1-exp(-(C*q+eps))) / (1-exp(-C))
    return min(bounds)