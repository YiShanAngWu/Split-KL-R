# 22 Nov, 2021
# Split kl for ternary random variables
# Consider Z\in\{z1,z2,z3\}, z_i\in\R
#
# Usage: python split-kl.py [RATE] [n]
#        RATE        : \in[0,1], p3 = RATE * (1-p2)
#        n           : number of samples    

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from math import log, sqrt, ceil
from tools import solve_kl_sup, splitkl, Emp_Bernstein, Unexp_Bernstein, Emp_MaV

PATH = 'plots/ternary/'
if not os.path.exists(PATH):
    os.makedirs(PATH)
RATE = float(sys.argv[1]) if len(sys.argv)>=2 else 0.5
n       = int(sys.argv[2]) if len(sys.argv)>=3 else 100

"""
Study the behavior of sum of n ternary random variables
taking values in {z1,z2,z3} with probability (p1,p2,p3).
"""
npoints = 1000 # number of different p2
window = 10 # report moving average

RV = np.array([-1, 0, 1]) # support of Z = {z1,z2,z3}
a, b = RV[0], RV[2] # range of Z
p2 = np.linspace(0., 1., npoints)
p3 = (1-p2) * RATE
p1 = (1-p2) * (1-RATE)
Probs = np.array([p1,p2,p3]).transpose() # probability distribution
delta = 0.05

""" compute the bounds """
kl_bound = np.zeros(npoints) # kl
EB_bound = np.zeros(npoints) # Empirical Bernstein
UB_bound = np.zeros(npoints) # Unexpected Bernstein
Skl_bound = np.zeros(npoints) # Split-kl


for i in range(npoints):
    # Define the probability distribution
    P = Probs[i]
    Ez = np.dot(P,RV) # true mean
    
    # Sample according to the probability distribution
    Sz = np.random.choice(RV, n, p=P)
    emp_Sz, var_Sz, z2 = Emp_MaV(Sz)
    
    # kl bound, rescale the r.v. to [0,1] to apply kl
    kl_bound[i] = (b-a)*solve_kl_sup((emp_Sz-a)/(b-a), log(1/delta)/n) + a - emp_Sz
    
    # Empirical Bernstein bound
    EB_bound[i] = Emp_Bernstein(emp_Sz, var_Sz, a, b, n, delta) - emp_Sz
    
    # Unexpected Bernstein bound
    UB_bound[i] = Unexp_Bernstein(emp_Sz, z2, b, n, delta) - emp_Sz

    # Split-kl
    mu = RV[1] # middle value
    Szp = np.maximum(np.zeros(n), Sz-mu)
    Szm = np.maximum(np.zeros(n), mu-Sz)
    Skl_bound[i] = splitkl(mu, a, b, np.mean(Szp), np.mean(Szm), log(2/delta)/n) - emp_Sz

""" Compute the moving average """
ma_kl = np.array([np.mean(kl_bound[i:i+window]) for i in range(npoints-window+1)])
ma_EB = np.array([np.mean(EB_bound[i:i+window]) for i in range(npoints-window+1)])
ma_UB = np.array([np.mean(UB_bound[i:i+window]) for i in range(npoints-window+1)])
ma_Skl = np.array([np.mean(Skl_bound[i:i+window]) for i in range(npoints-window+1)])

plt.rcParams.update({
    'font.size': 30,
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsfonts}'
})
""" Plot """
def plot():
    fig = plt.figure(figsize=(12,9.6))
    plt.plot(p2[:npoints-window+1], ma_kl, label='kl', linewidth=1.6)
    plt.plot(p2[:npoints-window+1], ma_EB, label='Emp-Bern', linewidth=1.6)
    plt.plot(p2[:npoints-window+1], ma_UB, label='Unexp-Bern', linewidth=1.6)
    plt.plot(p2[:npoints-window+1], ma_Skl, label='Split-kl', linewidth=1.6)
    title = r"$n=$"+str(n)+r", $p_1$="+str(RATE)+r"$(1-p_0)$"
    plt.title(title, fontsize=50)
    plt.ylabel(r"Bound $ - \hat{p}$", fontsize=50, fontweight='bold')
    plt.xlabel(r"$p_{0}$", fontsize=50, fontweight='bold')
    plt.xlim(0,1)
    plt.legend()
    plt.tight_layout()
    filename = "n"+str(n)+"_RATE"+str(RATE)
    path=PATH+filename+'.png'
    plt.savefig(path)
    plt.close()  
plot()
