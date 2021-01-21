import math 
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

def icl(n, alph):
  a_1 = stats.chi2.ppf(alph/2, df = n-1)
  a_2 = stats.chi2.ppf(1-alph/2, df = n-1)
  x = np.random.exponential(1,n)
  mx = np.mean(x)
  vx = np.var(x, ddof = 1)
  i_1, i_2 = (n-1)/a_2*vx , (n-1)/a_1*vx
  return mx, vx , (i_1 < 1 and i_2 >1), (i_2 - i_1), i_1, i_2

def interval(n,alph,m):
  obs = np.zeros((1000,6))
  count = 0
  for i in range(m):
    obs[i,:] = icl(n,alph)
  for i in range(m):
    if obs[i,2] == 1: 
      count +=1
  rho_hat = count/m
  return obs , rho_hat

def ic_tauxcouv(n,alph,m):
  a_1 = stats.t.ppf(alph/2,df = n-1)
  a_2 = stats.t.ppf(1-alph/2,df = n-1)
  rho_suite = interval(n,alph,m)[0][:,2]
  rho_hat = interval(n,alph,m)[1]
  rho_s = np.std(rho_suite, ddof = 1)
  i_1 = rho_hat - a_2*rho_s/math.sqrt(n-1)
  i_2 = rho_hat - a_1*rho_s/math.sqrt(n-1)
  return i_1, i_2

def ic_width(n,alph,m):
  a_1 = stats.t.ppf(alph/2,df = n-1)
  a_2 = stats.t.ppf(1-alph/2,df = n-1)
  w = interval(n,alph,m)[0][:,3]
  w_m = np.mean(w)
  w_s = np.std(w , ddof = 1)
  return (w_m - a_2*w_s/math.sqrt(n-1)) , (w_m - a_1*w_s/math.sqrt(n-1)), w_m

"""on a 1000 echantillons, chaque echantillon contient 100 observations."""

ic_tauxcouv(100,0.05,1000)

ic_width(100,0.05,1000)

q_1 = interval(100,0.05,1000) # la matrice contenant des data qu on va utiliser

q_1[1]

v = q_1[0][:,1] #obs

max(v) , min (v) #to know the length of axis

y_1 = q_1[0][:,4].tolist()
y_1.sort()
y_2 = q_1[0][:,5].tolist()
y_2.sort()
y_3 = q_1[0][:,1].tolist()
y_3.sort()

x = np.linspace(0,1,1000)
y_4 = stats.chi2.cdf(x , df = 99)
color = ["blue","green","orange","red"]
plt.plot(y_1,x)
plt.plot(y_2,x)
plt.plot(y_3,x)
plt.plot(chi2dist,x)
plt.legend(['IC left', 'IC right', 'Var estimator' , 'Chi-2'])
plt.grid()

"""1000 echantillon, chaque contient 1000 observations"""

ic_tauxcouv(1000,0.05,1000)

ic_width(1000,0.05,1000)

q_b = interval(1000,0.05,1000) # m = 1000

q_b[1]

v = q_b[0][:,1]
max(v) , min (v)

y_b1 = q_b[0][:,4].tolist()
y_b1.sort()
y_b2 = q_b[0][:,5].tolist()
y_b2.sort()
y_b3 = q_b[0][:,1].tolist()
y_b3.sort()

x = np.linspace(0,1,1000)
y_b4 = stats.chi2.ppf(x , df = 999)/999
color = ["blue","green","orange","red"]
plt.plot(y_b1,x)
plt.plot(y_b2,x)
plt.plot(y_b3,x)
plt.plot(y_b4,x)
plt.title("with n = 1000")
plt.legend(['IC left', 'IC right', 'Var estimator' , 'Chi-2'])
plt.grid()