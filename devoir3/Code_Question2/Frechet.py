# -*- coding: utf-8 -*-
"""
Question 2 b) et c)
Pour la méthode de Fréchet
"""
import math
import numpy as np
from numpy.linalg import cholesky
import scipy.stats as stats
import matplotlib.pyplot as plt

def frechet_1(rho_s):# this one shows the related parameter value
    p=(rho_s+1)/2
    print("parameter p of the frechet copula=",p)
    p_test=np.random.uniform(0,1)
    u_1=np.random.uniform(0,1)
    u_2=0
    if p_test<p:
        u_2=u_1
    else:
        u_2=1-u_1
    return u_1,u_2

def frechet(rho_s):# the one used in simulation
    p=(rho_s+1)/2
    #print("parameter p of the frechet copula=",p)
    p_test=np.random.uniform(0,1)
    u_1=np.random.uniform(0,1)
    u_2=0
    if p_test<p:
        u_2=u_1
    else:
        u_2=1-u_1
    return u_1,u_2

#question b
j=0
for i in range(1000000):
    U=frechet(0.9)
    if U[0]+U[1]<0.1:#U is a tuple here
        j+=1
j/1e6
#for N=1e3, j=0.047
#for N=1e4, j=0.0494
#for N=1e5, j=0.047
#for N=1e6, j=0.047785

#question c
j=0
for i in range(1000000):
    U=frechet(0.8)
    if U[0]-U[1]<0.05:#U is a tuple here
        j+=1
j/1e6
#for N=1e3, j=0.945
#for N=1e4, j=0.9551
#for N=1e5, j=0.95333
#for N=1e6, j=0.952566




