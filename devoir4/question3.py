# -*- coding: utf-8 -*-
"""
Created on Wed Nov 9 12:47:05 2020
"""
import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

###question b&d
def desite(x,y):
    re = -x*x-y*y
    f = math.exp(re/2)/(2*math.pi)
    return f

def generator (a,b,c):
    u_1,u_2,u_3 = np.random.uniform(), np.random.uniform(), np.random.uniform()
    x = a*u_1
    y = b*u_2
    z = c*u_3
    while desite(x,y) <= z : #accept le premier Z pour Z <= f(x,y) 
        u_1,u_2,u_3 = np.random.uniform(), np.random.uniform(), np.random.uniform()
        x = a*u_1
        y = b*u_2
        z = c*u_3
    return x,y

def simuler(a,b,n):
    i = 0
    m = np.zeros((n,2))
    while i < n:
        m[i][0], m[i][1] = generator(a,b,1/(math.pi*2))
        i+= 1
    return m

m = simuler(4,5,1000)
plt.scatter(m[:,0],m[:,1],s=12)
plt.xlim([0,4])
plt.ylim([0,5]) 
plt.xlabel('X')
plt.ylabel('Y')           
plt.grid()  

      
###question e
stats.norm.ppf(0.95)

def normgene (a,b):
    u_1,u_2 = np.random.uniform(), np.random.uniform() 
    x = stats.norm.ppf(u_1)
    y = stats.norm.ppf(u_2)
    while x<0 or x>a:
        x = stats.norm.ppf(np.random.uniform())
    while y<0 or y>b:
        y = stats.norm.ppf(np.random.uniform())
    return x,y

def simuler_2(a,b,n):
    i = 0
    m = np.zeros((n,2))
    while i < n:
        m[i][0], m[i][1] = normgene(a,b)
        i +=1
    return m
    
m_2 = simuler_2(4,5,1000)
plt.scatter(m_2[:,0],m_2[:,1],s=10)       
                
plt.subplot(1,2,1)  
plt.xlim([0,4])
plt.ylim([0,5]) 
plt.scatter(m[:,0],m[:,1],s=12)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('methode 1')
plt.grid()
plt.subplot(1,2,2)
plt.xlim([0,4])
plt.ylim([0,5]) 
plt.scatter(m_2[:,0],m_2[:,1],s=10)  
plt.xlabel('X')
plt.ylabel('Y')           
plt.title('methode 2')
plt.grid()        