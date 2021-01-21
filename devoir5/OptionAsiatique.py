from google.colab import files
uploaded = files.upload()

"""the standard uniform random numbers we used in this exemple is generated from the generator MRG32k3a"""

import numpy as np
import scipy.stats as stats
import math

randomnum = np.loadtxt('randomdouble_MRG32k3a.txt')

"""import the fixed parameters"""

sigma = 0.2 #all the fixed parameters
r = 0.08
s_0 = 100
T = 120/365
n = 10000
K = np.array([80,90,100,110])

"""creat the arrays to contain the time

time case 1
"""

t_1 = np.ones(10) #time (1)
for j in range(10):
  t_1[j] = (110 + j+1)/365

"""time case 2"""

t_2 = np.ones(10) #time (2)
for j in range(10):
  t_2[j] = (12*(j+1))

"""time case 3"""

t_3 = np.ones(120) # time (3)
for j in range (120):
  t_3[j] = (j+1)/365

def price(t, d) : # prepare this for methode 5
  b = np.zeros(d)
  s = np.ones(d)
  b2 = np.zeros(d)
  s2 = np.ones(d)
  for i in range(1,d):
    u = np.random.uniform(0,1)
    b[i] = b[i-1] + math.sqrt(t[i]-t[i-1])*stats.norm.ppf(u)
    s[i-1] = s_0*math.exp((r-sigma**2/2)*t[i] + sigma* b[i-1])
    b2[i] = b2[i-1] + math.sqrt(t[i]-t[i-1])*stats.norm.ppf(1-u)
    s2[i] = s_0*math.exp((r-sigma**2/2)*t[i] + sigma* b2[i])
  s[d-1] = s_0*math.exp((r-sigma**2/2)*t[d-1] + sigma* b[d-1])
  s2[d-1] = s_0*math.exp((r-sigma**2/2)*t[d-1] + sigma* b2[d-1])
  somme_0 = np.sum(s)
  multi_0 = np.prod(s)
  sum_anti = np.sum(s2)
  multi_anti = np.prod(s2)
  return somme_0, multi_0, sum_anti, multi_anti

"""the following fonction is for question 1 to 4!!"""

def for124(t, d) :
  b = np.zeros(d)
  s = np.ones(d)
  for i in range(1,d):
    b[i] = b[i-1] + math.sqrt(t[i]-t[i-1])*stats.norm.ppf(randomnum[i])
    s[i-1] = s_0*math.exp((r-sigma**2/2)*t[i] + sigma* b[i-1])
  s[d-1] = s_0*math.exp((r-sigma**2/2)*t[d-1] + sigma* b[d-1])
  somme = np.sum(s)
  multi = np.prod(s)
  return somme, multi

def for124(t, d) :
  b = np.zeros(d)
  s = np.ones(d)
  for i in range(1,d):
    u = np.random.uniform(0,1)
    b[i] = b[i-1] + math.sqrt(t[i]-t[i-1])*stats.norm.ppf(u)
    s[i-1] = s_0*math.exp((r-sigma**2/2)*t[i] + sigma* b[i-1])
  s[d-1] = s_0*math.exp((r-sigma**2/2)*t[d-1] + sigma* b[d-1])
  somme = np.sum(s)
  multi = np.prod(s)
  return somme, multi

def cal(t,d,n,k):
  arith = np.ones(n)
  geo = np.ones(n)
  sum = np.ones(n)
  for i in range (n):
    re = for124(t,d)
    sum[i] = re[0]
    arith[i] = max(0, math.exp(-r*T))*(re[0]/d - k)
    geo[i] = max(0, math.exp(-r*T))*(re[1]**(1/d)-k)
  m_a = np.mean(arith)
  return arith, geo, m_a, sum

"""Calculate the CV1"""

def cv1(t,d,k):
  timecv1 = np.sum(t)
  timecv2 = 0
  for i in range(1,d):
    timecv2 += (t[i]-t[i-1])*(d-i+1)**2
  my = math.log(s_0) + (r - math.pow(sigma,2)/2)/d*timecv1
  sy2 = math.pow((sigma/d),2)*timecv2
  d = (-math.log(k)+my)/sy2
  cv1 = math.exp(-r*T)*(math.exp(my + sy2/2)*stats.norm.cdf(d+math.sqrt(sy2))-k*stats.norm.cdf(d))
  return cv1

"""Calculate the CV2"""

def cv2(t,d,k):
  timecv2 = 0
  for i in range(d):
    timecv2 += math.exp(r*t[i])
  cv2 = s_0*timecv2
  return cv2

"""calculate all the variance/covariance"""

def calv(t,d,n,k):
  arith, geo, mA, sum = cal(t,d,n,k)
  controlv_1 = cv1(t,d,k)
  controlv_2 = cv2(t,d,k)
  varA = varG = varsum =covAG = covAsum = covGsum = 0
  for i in range(n):
    varA += math.pow((arith[i]-mA),2)/(n-1)
    varG += math.pow((geo[i]-controlv_1),2)/(n-1)
    varsum += math.pow((sum[i] -controlv_2),2)/(n-1)
    covAG += (geo[i]-controlv_1)*(arith[i]-mA)/(n-1)
    covAsum += (sum[i]-controlv_2)*(arith[i]-mA)/(n-1)
    covGsum += (sum[i]-controlv_2)*(geo[i]-controlv_1)/(n-1)
  return varA, varG, varsum, covAG, covAsum, covGsum, arith, geo, controlv_1, controlv_2, sum

#calv(t_1,10,10,80)

"""Methode 1: Withour any VRT"""

def noVRT(t,d,n,k):
  arith = np.ones(n)
  for i in range (n):
    re = for124(t,d)
    arith[i] = max(0, math.exp(-r*T))*(re[0]/d - k)
  return arith

"""Methode 2 : Geometric Variance Control"""

def vcGeo(t,d,n,k):
  v = calv(t,d,n,k)
  varA, varG, varsum, covAG, covAsum, covGsum, arith, geo, controlv_1, controlv_2, sum = v
  betaG = covAG/varG
  re_geo = np.ones(n)
  for i in range(n):
    re_geo[i] = arith[i]-betaG*(geo[i]-controlv_1)
  return re_geo

"""Methode 3 : Sum Control Variable"""

def vcSum(t,d,n,k):
  v = calv(t,d,n,k)
  varA, varG, varsum, covAG, covAsum, covGsum, arith, geo, controlv_1, controlv_2, sum = v
  beta = covAsum/varsum
  re_sum = np.ones(n)
  for i in range(n):
    re_sum[i] = arith[i]-beta*(sum[i]-controlv_2)
  return re_sum

"""Methode 4 : Sum of two control variable"""

def twoVC(t,d,n,k):
  v = calv(t,d,n,k)
  varA, varG, varsum, covAG, covAsum, covGsum, arith, geo, controlv_1, controlv_2, sum = v
  beta_1 = (varsum*covAG - covGsum*covAsum)/(varG*varsum - covGsum**2)
  beta_2 = (varG*covAsum - covGsum*covAG)/(varG*varsum - covGsum**2)
  re_2vc = np.ones(n)
  for i in range(n):
    re_2vc[i] = arith[i]-beta_1*(geo[i]-controlv_1)-beta_2*(sum[i]-controlv_2)
  return re_2vc

"""Methode 5 : with both control variables and antithetic variates"""

def methode5(t, d) :
  b = np.zeros(d)
  s = np.ones(d)
  b2 = np.zeros(d)
  s2 = np.ones(d)
  for i in range(1,d):
    u = np.random.uniform(0,1)
    b[i] = b[i-1] + math.sqrt(t[i]-t[i-1])*stats.norm.ppf(u)
    s[i-1] = s_0*math.exp((r-sigma**2/2)*t[i] + sigma* b[i-1])
    b2[i] = b2[i-1] + math.sqrt(t[i]-t[i-1])*stats.norm.ppf(1-u)
    s2[i] = s_0*math.exp((r-sigma**2/2)*t[i] + sigma* b2[i])
  s[d-1] = s_0*math.exp((r-sigma**2/2)*t[d-1] + sigma* b[d-1])
  s2[d-1] = s_0*math.exp((r-sigma**2/2)*t[d-1] + sigma* b2[d-1])
  somme_0 = np.sum(s)
  multi_0 = np.prod(s)
  sum_anti = np.sum(s2)
  multi_anti = np.prod(s2)
  return somme_0, multi_0, sum_anti, multi_anti

def cal5(t,d,n,k):
  arith = np.ones(n)
  geo = np.ones(n)
  sum = np.ones(n)
  for i in range (n):
    re = methode5(t,d)
    somme_0, multi_0, sum_anti, multi_anti = re
    arith[i] = (max(0, math.exp(-r*T)*(somme_0/d - k))+max(0,math.exp(-r*T)*(sum_anti/d - k)))/2
    geo[i] = (max(0, math.exp(-r*T)*(multi_0**(1/d)-k))+max(0,math.exp(-r*T)*(multi_anti**(1/d) - k)))/2
    sum[i] = (somme_0 + sum_anti)/2
  ma = np.mean(arith)
  return arith, geo, ma, sum

def par5(t,d,n,k):
  arith, geo, mA, sum = cal5(t,d,n,k)
  controlv_1 = cv1(t,d,k)
  controlv_2 = cv2(t,d,k)
  varA = varG = varsum =covAG = covAsum = covGsum = 0
  for i in range(n):
    varA += math.pow((arith[i]-mA),2)/(n-1)
    varG += math.pow((geo[i]-controlv_1),2)/(n-1)
    varsum += math.pow((sum[i] -controlv_2),2)/(n-1)
    covAG += (geo[i]-controlv_1)*(arith[i]-mA)/(n-1)
    covAsum += (sum[i]-controlv_2)*(arith[i]-mA)/(n-1)
    covGsum += (sum[i]-controlv_2)*(geo[i]-controlv_1)/(n-1)
  return varA, varG, varsum, covAG, covAsum, covGsum, arith, geo, controlv_1, controlv_2, sum

def anti2vc(t,d,n,k):
  v = par5(t,d,n,k)
  varA, varG, varsum, covAG, covAsum, covGsum, arith, geo, controlv_1, controlv_2, sum = v
  beta_1 = (varsum*covAG - covGsum*covAsum)/(varG*varsum - covGsum**2)
  beta_2 = (varG*covAsum - covGsum*covAG)/(varG*varsum - covGsum**2)
  re_anti2vc = np.ones(n)
  for i in range(n):
    re_anti2vc[i] = arith[i]-beta_1*(geo[i]-controlv_1)-beta_2*(sum[i]-controlv_2)
  return re_anti2vc

""" time (1) avec k =80"""

re_novrt_1 = noVRT(t_1,10,10000,80)

vrfbase_11 = np.var(re_novrt_1)  #VRF base!

vrfbase_11

re_geo_1 = vcGeo(t_1,10,10000,80)

vrfbase_11/np.var(re_geo_1)   ##vrf geo 11

re_sum_1 = vcSum(t_1,10,10000,80)

re_2vc_1 = twoVC(t_1,10,10000,80)

re_anti_11 = anti2vc(t_1,10,10000,80) #methode 5

vrfbase_11/np.var(re_sum_1)  ##vrf Sum 11

vrfbase_11/np.var(re_2vc_1) # vrf 2VC 11

vrfbase_11/np.var(re_anti_11) # vrf method5

a,b,c,d,e = np.mean(re_novrt_1), np.mean(re_geo_1), np.mean(re_sum_1), np.mean(re_2vc_1), np.mean(re_anti_11)

a,b,c,d,e

"""time(1) avec k =90"""

re_novrt_2 = noVRT(t_1,10,10000,90)
re_geo_2 = vcGeo(t_1,10,10000,90)
re_sum_2 = vcSum(t_1,10,10000,90)
re_2vc_2 = twoVC(t_1,10,10000,90)

re_anti_12 = anti2vc(t_1,10,10000,90) #methode 5

vrfbase_12 = np.var(re_novrt_2)  #VRF base!

vrfbase_12

vrfbase_12/np.var(re_geo_2)   ##vrf geo 12

vrfbase_12/np.var(re_sum_2)   ##vrf sum 12

vrfbase_12/np.var(re_2vc_2)   ##vrf both vc 12

vrfbase_12/np.var(re_anti_12) # vrf method5 t=1 k=90

a1,b1,c1,d1,e1 = np.mean(re_novrt_2), np.mean(re_geo_2), np.mean(re_sum_2), np.mean(re_2vc_2), np.mean(re_anti_12) ###toooooo run

a1,b1,c1,d1,e1

"""time(1) avec k =100"""

re_novrt_3 = noVRT(t_1,10,10000,100)
re_geo_3 = vcGeo(t_1,10,10000,100)
re_sum_3 = vcSum(t_1,10,10000,100)
re_2vc_3 = twoVC(t_1,10,10000,100)

re_anti_13 = anti2vc(t_1,10,10000,100) #methode 5

vrfbase_13 = np.var(re_novrt_3)  #VRF base!

vrfbase_13

vrfbase_13/np.var(re_geo_3)   ##vrf geo 13

vrfbase_13/np.var(re_sum_3)   ##vrf sum 13

vrfbase_13/np.var(re_2vc_3)   ##vrf both vc 13

vrfbase_13/np.var(re_anti_13) # vrf method5 t=1 k=100

a2,b2,c2,d2,e2 = np.mean(re_novrt_3), np.mean(re_geo_3), np.mean(re_sum_3), np.mean(re_2vc_3), np.mean(re_anti_13) ###toooooo run

a2,b2,c2,d2,e2

"""time(1) avec k =110"""

re_novrt_4 = noVRT(t_1,10,10000,110)
re_geo_4 = vcGeo(t_1,10,10000,110)
re_sum_4 = vcSum(t_1,10,10000,110)
re_2vc_4 = twoVC(t_1,10,10000,110)

re_anti_14 = anti2vc(t_1,10,10000,110) #methode 5

vrfbase_14 = np.var(re_novrt_4)  #VRF base!

vrfbase_14

vrfbase_14/np.var(re_geo_4)   ##vrf geo 14

vrfbase_14/np.var(re_sum_4)   ##vrf sum 14

vrfbase_14/np.var(re_2vc_4)   ##vrf both vc 14

vrfbase_14/np.var(re_anti_14) # vrf method5 t=1 k=110

a4,b4,c4,d4 = np.mean(re_novrt_4), np.mean(re_geo_4), np.mean(re_sum_4), np.mean(re_2vc_4) ###toooooo run

a4,b4,c4,d4

"""time(2) avec k = 80"""

re_novrt_21 = noVRT(t_2,10,10000,80)
re_geo_21 = vcGeo(t_2,10,10000,80)
re_sum_21 = vcSum(t_2,10,10000,80)
re_2vc_21 = twoVC(t_2,10,10000,80)
re_anti_21 = anti2vc(t_2,10,10000,80)

re_novrt_21 ##### not gonna use this just for checking

vrfbase_21 = np.var(re_novrt_21)  #VRF base!

vrfbase_21

vrfbase_21/np.var(re_geo_21)   ##vrf geo 21 #t = 2, k = 80

vrfbase_21/np.var(re_sum_21)   ##vrf sum 21

vrfbase_21/np.var(re_2vc_21)   ##vrf both vc 21

vrfbase_21/np.var(re_anti_21) # vrf method5 t=2 k=80

a,b,c,d,e = np.mean(re_novrt_21), np.mean(re_geo_21), np.mean(re_sum_21), np.mean(re_2vc_21), np.mean(re_anti_21) ###toooooo run

a,b,c,d,e

"""time(2) avec k = 90"""

re_novrt_22 = noVRT(t_2,10,10000,90)
re_geo_22 = vcGeo(t_2,10,10000,90)
re_sum_22 = vcSum(t_2,10,10000,90)
re_2vc_22 = twoVC(t_2,10,10000,90)
re_anti_22 = anti2vc(t_2,10,10000,90) #methode 5

vrfbase_22 = np.var(re_novrt_21)  #VRF base!

vrfbase_22

vrfbase_22/np.var(re_geo_22)   ##vrf geo 22 # t = ii , k = 90

vrfbase_22/np.var(re_sum_22)   ##vrf sum 22

vrfbase_22/np.var(re_2vc_22)   ##vrf both vc 22

vrfbase_22/np.var(re_anti_22) # vrf method5 t=2 k=90

a,b,c,d,e = np.mean(re_novrt_22), np.mean(re_geo_22), np.mean(re_sum_22), np.mean(re_2vc_22), np.mean(re_anti_22) ###toooooo run

a,b,c,d,e

"""time(2) avec k = 100"""

re_novrt_23 = noVRT(t_2,10,10000,100)
re_geo_23 = vcGeo(t_2,10,10000,100)
re_sum_23 = vcSum(t_2,10,10000,100)
re_2vc_23 = twoVC(t_2,10,10000,100)
re_anti_23 = anti2vc(t_2,10,10000,100)

vrfbase_23 = np.var(re_novrt_23)  #VRF base!

vrfbase_23

vrfbase_23/np.var(re_geo_23)   ##vrf geo 23  ### t =

vrfbase_23/np.var(re_sum_23)   ##vrf sum 23

vrfbase_23/np.var(re_2vc_23)   ##vrf both vc 23

vrfbase_23/np.var(re_anti_23) # vrf method5 t=2 k=100

a,b,c,d,e = np.mean(re_novrt_23), np.mean(re_geo_23), np.mean(re_sum_23), np.mean(re_2vc_23), np.mean(re_anti_23) ###toooooo run

a,b,c,d,e

"""time(2) avec k = 110"""

re_novrt_24 = noVRT(t_2,10,10000,110)
re_geo_24 = vcGeo(t_2,10,10000,110)
re_sum_24 = vcSum(t_2,10,10000,110)
re_2vc_24 = twoVC(t_2,10,10000,110)

re_anti_24 = anti2vc(t_2,10,10000,110)

vrfbase_24 = np.var(re_novrt_24)  #VRF base!

vrfbase_24

vrfbase_24/np.var(re_geo_24)   ##vrf geo 24

vrfbase_24/np.var(re_sum_24)   ##vrf sum 24

vrfbase_24/np.var(re_2vc_24)   ##vrf both vc 24

vrfbase_24/np.var(re_anti_24) # vrf method5 t=2 k=110

a,b,c,d,e = np.mean(re_novrt_24), np.mean(re_geo_24), np.mean(re_sum_24), np.mean(re_2vc_24), np.mean(re_anti_24) ###toooooo run

a,b,c,d,e

"""time(3) avec k = 80"""

re_novrt_31 = noVRT(t_3,120,10000,80)
re_geo_31 = vcGeo(t_3,120,10000,80)
re_sum_31 = vcSum(t_3,120,10000,80)
re_2vc_31 = twoVC(t_3,120,10000,80)

re_anti_31 = anti2vc(t_3,120,10000,80)    ####Running!!!!!!!!!!!!!!!!!!

vrfbase_31 = np.var(re_novrt_31)  #VRF base!

vrfbase_31

vrfbase_31/np.var(re_geo_31)   ##vrf geo 31

vrfbase_31/np.var(re_sum_31)   ##vrf sum 31

vrfbase_31/np.var(re_2vc_31)   ##vrf both vc 31

vrfbase_31/np.var(re_anti_31) # vrf method5 t=3 k=80

a,b,c,d,e = np.mean(re_novrt_31), np.mean(re_geo_31), np.mean(re_sum_31), np.mean(re_2vc_31), np.mean(re_anti_31) ###toooooo run

a,b,c,d,e  #### to edit

"""time (3) avec k = 90"""

re_novrt_32 = noVRT(t_3,120,10000,90)
re_geo_32 = vcGeo(t_3,120,10000,90)
re_sum_32 = vcSum(t_3,120,10000,90)
re_2vc_32 = twoVC(t_3,120,10000,90)

re_anti_32 = anti2vc(t_3,120,10000,90)

vrfbase_32 = np.var(re_novrt_32)  #VRF base!

vrfbase_32

vrfbase_32/np.var(re_geo_32)   ##vrf geo 32

vrfbase_32/np.var(re_sum_32)   ##vrf sum 32

vrfbase_32/np.var(re_2vc_32)   ##vrf both vc 32

vrfbase_32/np.var(re_anti_32) # vrf method5 t=3 k=90

a,b,c,d,e = np.mean(re_novrt_32), np.mean(re_geo_32), np.mean(re_sum_32), np.mean(re_2vc_32), np.mean(re_anti_32) ###toooooo run

a,b,c,d,e

"""time (3) avec k = 100"""

re_novrt_33 = noVRT(t_3,120,10000,100)
re_geo_33 = vcGeo(t_3,120,10000,100)
re_sum_33 = vcSum(t_3,120,10000,100)
re_2vc_33 = twoVC(t_3,120,10000,100)

re_anti_33 = anti2vc(t_3,120,10000,100)

vrfbase_33 = np.var(re_novrt_33)  #VRF base!

vrfbase_33

vrfbase_33/np.var(re_geo_33)   ##vrf geo 33

vrfbase_33/np.var(re_sum_33)   ##vrf sum 33

vrfbase_33/np.var(re_2vc_33)   ##vrf both vc 33

vrfbase_33/np.var(re_anti_33) # vrf method5 t=3 k=100

a,b,c,d,e = np.mean(re_novrt_33), np.mean(re_geo_33), np.mean(re_sum_33), np.mean(re_2vc_33), np.mean(re_anti_33) ###toooooo run

a,b,c,d,e

"""time (3) avec k = 110"""

re_novrt_34 = noVRT(t_3,120,10000,110)
re_geo_34 = vcGeo(t_3,120,10000,110)
re_sum_34 = vcSum(t_3,120,10000,110)
re_2vc_34 = twoVC(t_3,120,10000,110)

re_anti_34 = anti2vc(t_3,120,10000,110)

vrfbase_34 = np.var(re_novrt_34)  #VRF base!

vrfbase_34

vrfbase_34/np.var(re_geo_34)   ##vrf geo 34

vrfbase_34/np.var(re_sum_34)   ##vrf sum 34

vrfbase_34/np.var(re_2vc_34)   ##vrf both vc 34

vrfbase_34/np.var(re_anti_34) # vrf method5 t=3 k=110

a,b,c,d,e = np.mean(re_novrt_34), np.mean(re_geo_34), np.mean(re_sum_34), np.mean(re_2vc_34), np.mean(re_anti_34) ###toooooo run

a,b,c,d,e