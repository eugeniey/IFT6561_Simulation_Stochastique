# -*- coding: utf-8 -*-
"""
Question 2 b) 
Pour la méthode de Norta
"""
import math
import numpy as np
from numpy.linalg import cholesky
import scipy.stats as stats
import matplotlib.pyplot as plt

def norta(rho_s):
    rho=2*np.sin(math.pi*rho_s/6)
    sigma=np.mat([[1,rho],[rho,1]])
    L=cholesky(sigma)
    Z=np.mat([np.random.normal(size=100),np.random.normal(size=100)])
    Y=np.dot(L,Z)
    U_1,U_2=stats.norm.cdf(Y[0]),stats.norm.cdf(Y[1])
    return U_1,U_2

def norta_b(rho_s):
    rho=2*np.sin(math.pi*rho_s/6)
    sigma=np.mat([[1,rho],[rho,1]])
    L=cholesky(sigma)
    Z=np.mat([np.random.normal(size=1),np.random.normal(size=1)])
    Y=np.dot(L,Z)
    U_1,U_2=stats.norm.cdf(Y[0]),stats.norm.cdf(Y[1])
    return U_1,U_2
#plt.hist(U_1) if want to verify the distribution of U
    
#question b
j=0
for i in range(1000000):
    U=norta_b(0.9)
    if U[0][0][0]+U[1][0][0]<0.1:#U is a tuple here
        j+=1
j/1e6
#result j=42673
#P=j/1e6=0.042673

#question c
j=0
for i in range(1000000):
    U=norta_b(0.8)
    if U[0][0][0]-U[1][0][0]<0.05:
        j+=1
j/1e6
#resultat j=639009
#P=j/1e6=0.639009
#j'ai essayé avec différentes valeurs de N(a part de 1e6)
#le plus grand quand N augmente
#le plus grand la probabilité devient

        