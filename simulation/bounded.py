# 19 May, 2022
# Split-kl for random variables in [0,1] using beta distribution with parameters alpha, beta.
#
# Usage: python split-kl-general.py [SETTING] [n]
#        Setting     : 'ConstM', 'Spectrum'
#                      'ConstM'        = distributions with const mean p and different variances
#                      'Spectrum'      = distributions with mean p taken in [0,1]
#        n           : number of samples    

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from math import log, sqrt, ceil
from tools import solve_kl_sup, splitkl, Emp_Bernstein, Unexp_Bernstein, Emp_MaV

PATH = 'plots/bounded/'
if not os.path.exists(PATH):
    os.makedirs(PATH)

SETTING = sys.argv[1] if len(sys.argv)>=2 else 'ConstM'
n       = int(sys.argv[2]) if len(sys.argv)>=3 else 100

"""
Study the behavior of sum of n bounded random variables
taking values in [a,b].
Consider a beta distribution with parameters alpha, beta.
"""
npoints = 1000 # the number of distributions to test
window = 10 # report moving average

if SETTING == 'Spectrum':
    """ Consider the combination of two spectrums:
        1. \beta=5, \alpha\in[0.01,5]
        2. \alpha=5, \beta\in[0.01,5]
        The overall pdf is a bell-shape, moving from left to right.
    """
    Prob1 = np.array([[0.01*(i+1),5] for i in range(int(npoints/2))])
    Prob2 = np.array([[5,5-0.01*i] for i in range(int(npoints/2))])
    Prob = np.concatenate((Prob1,Prob2), axis=0)
elif SETTING == 'ConstM':
    """ Consider distributions with same mean but different variance:
        alpha=beta\in[0.01,10]
    """
    Prob = np.array([[0.01*(i+1),0.01*(i+1)] for i in range(npoints)])
else:
    print("No such SETTING!")
    sys.exit(1)

P = Prob[:,0]/(Prob[:,0]+Prob[:,1]) # true mean of beta distribution
Var = (Prob[:,0]*Prob[:,1])/((Prob[:,0]+Prob[:,1])**2*(Prob[:,0]+Prob[:,1]+1))
a = 0
b = 1
delta = 0.05

""" compute the bounds """
kl_bound = np.zeros(npoints) # kl
EB_bound = np.zeros(npoints) # Empirical Bernstein
UB_bound = np.zeros(npoints) # Unexpected Bernstein
Skl_bound = np.zeros(npoints) # Split-kl


for i in range(npoints):
    # Define the probability distribution
    Sz = np.random.beta(a=Prob[i,0], b=Prob[i,1], size=n)
    emp_Sz, var_Sz, z2 = Emp_MaV(Sz)
    
    # kl bound, rescale the r.v. to [0,1] to apply kl
    kl_bound[i] = (b-a)*solve_kl_sup((emp_Sz-a)/(b-a), log(1/delta)/n) + a - emp_Sz
    
    # Empirical Bernstein bound
    EB_bound[i] = Emp_Bernstein(emp_Sz, var_Sz, a, b, n, delta) - emp_Sz
    
    # Unexpected Bernstein bound
    UB_bound[i] = Unexp_Bernstein(emp_Sz, z2, b, n, delta) - emp_Sz

    # Split-kl
    mu = 0.5 # middle value
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

def plot():
    fig = plt.figure(figsize=(12,9.6))
    if SETTING == 'Spectrum':
        plt.plot(P[:npoints-window+1], ma_kl, label='kl', linewidth=1.6)
        plt.plot(P[:npoints-window+1], ma_EB, label='Emp-Bern', linewidth=1.6)
        plt.plot(P[:npoints-window+1], ma_UB, label='Unexp-Bern', linewidth=1.6)
        plt.plot(P[:npoints-window+1], ma_Skl, label='Split-kl', linewidth=1.6)
        plt.xlabel(r"$p$", fontsize=50, fontweight='bold')
    elif SETTING == 'ConstM':
        plt.plot(Var[:npoints-window+1], ma_kl, label='kl', linewidth=1.6)
        plt.plot(Var[:npoints-window+1], ma_EB, label='Emp-Bern', linewidth=1.6)
        plt.plot(Var[:npoints-window+1], ma_UB, label='Unexp-Bern', linewidth=1.6)
        plt.plot(Var[:npoints-window+1], ma_Skl, label='Split-kl', linewidth=1.6)
        plt.xlabel(r"$\mathbb{V}[Z]$", fontsize=50, fontweight='bold')
    
    title = r"$n="+str(n)+r"$"
    plt.title(title, fontsize=50)
    plt.ylabel(r"Bound $ - \hat{p}$", fontsize=50, fontweight='bold')
    plt.legend()
    plt.tight_layout()
    filename = "n"+str(n)+"_"+SETTING
    path=PATH+filename+'.png'
    plt.savefig(path)
    plt.close()  
plot()
