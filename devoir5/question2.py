import math 
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

x = np.random.exponential(1,100)

"""La fct boots() est pour rééchantilloner les X_n ainsi obtenir les X_n*"""

def boots(x,m):
  y = np.zeros((100,m))
  for i in range(m):
    for j in range (100):
      b = np.random.randint(0,99)
      y[j,i] = x[b]
  return y

"""we want the Y_n from X_n*, then we put them in order(small to large)"""

def Kn(x_star):
  y = np.zeros(np.size(x_star,1))
  for i in range (np.size(x_star,1)):
    y[i] = np.var(x_star[:,i], ddof = 1)
  y.sort()
  return y

"""we know the theoratical variance of X_n = 1, we use this value to count the number of effective confidence interval"""

def ic (m,alp,n):
  num_1 = int(m*(1-(alp/2)))
  num_2 = int(m*(alp/2))
  i_inf, i_sup = np.zeros(n), np.zeros(n)
  sn, w = np.zeros(n), np.zeros(n)
  count = 0
  rho_data = np.zeros(n)
  for j in range(n):
    x = np.random.exponential(1,100)
    x_star = boots(x,m)
    yn = Kn(x_star)
    sn[j] = np.var(x, ddof =1)
    i_inf[j] = 2*sn[j] - yn[num_1]
    i_sup[j] = 2*sn[j] - yn[num_2]
    w[j] = i_sup[j] - i_inf[j]
  for i in range(n):
    if i_inf[i]<=1 and i_sup[i]>=1:
      count +=1
      rho_data[i] = 1
  e_w = np.sum(w)/n
  print("rho =", count/n)
  print("espérance W = ", e_w)    
  return i_inf, i_sup, w, sn, (count/n), rho_data

def ic_tauxcouv(m,alp,n): #等下删
  a_1 = stats.t.ppf(alp/2,df = n-1)
  a_2 = stats.t.ppf(1-alp/2,df = n-1)
  rho_suite = ic(m,alp,n)[5]
  rho_hat = ic(m,alp,n)[4]
  rho_s = np.std(rho_suite, ddof = 1)
  i_1 = rho_hat - a_2*rho_s/math.sqrt(n-1)
  i_2 = rho_hat - a_1*rho_s/math.sqrt(n-1)
  return i_1, i_2

re = ic(1000,0.05,1000)

ic_tauxcouv(100,0.05,1000)

"""Interpreter en graphique"""

x1 = np.linspace(0,1,1000)
y1 = np.sort(re[0])
y2 = np.sort(re[1])
y3 = np.sort(re[3])
color = ["blue","green","orange"]
plt.plot(y1,x1)
plt.plot(y2,x1)
plt.plot(y3,x1)
plt.legend(['IC left', 'IC right', 'Var estimator' ])
plt.grid()