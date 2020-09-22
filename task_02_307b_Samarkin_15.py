from scipy.constants import c, pi
from scipy.special import spherical_jn as jn, spherical_yn as yn
import numpy as np
from pathlib import Path as pth
import matplotlib.pyplot as plt

#Тестовые переменные
D = 15e-3
fmin = 0.01e9
fmax = 45e9

f = np.linspace(fmin, fmax, 400)
r = D/2

def hn(n, x):
    return complex(jn(n, x), yn(n,x))

def bn(n, x):
    upr_n = x * jn(n-1, x) - n*jn(n, x)
    lwr_n = x * hn(n-1, x) - n * hn(n, x)
    return upr_n / lwr_n

def an(n, x):
    return jn(n, x) / hn(n, x) 

def f_sigma(x):
    lmbd = c / x #Лямбда
    k = 2 * pi / lmbd #Волновое число
    kr = k*r
    Sum = 0+0j
    for n in range(1, 45):
        Sum += (-1 ** n) * (n + 1/2) * (bn(n, kr) - an(n, kr))
    return ((lmbd**2) / pi) * (np.abs(Sum)**2) 

sigma = [f_sigma(x) for x in f]

p = pth('results')
res = p / 'task_02_307b_Samarkin_15.txt'
if not p.exists():
  p.mkdir(exist_ok = True)

if p.exists():
  with res.open('w') as fl:
    i = 0
    for x in f:
      fl.write('%11.10f    %.10f\n' % (x, sigma[i]))
      i += 1
else:
  print("Error: \'results\' dir does not exists - access denied")


plt.plot(f/1e9, sigma)
plt.xlabel('f, ГГц')
plt.ylabel('$\sigma, м^2$')
plt.grid()
plt.show()

