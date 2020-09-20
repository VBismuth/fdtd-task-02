from scipy.constants import c, pi
from scipy.special import spherical_jn as jn, spherical_yn as yn
import numpy as np
import pathlib as pth
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
    if lwr_n != 0: return upr_n / lwr_n
    else: return 0

def an(n, x):
    if hn(n, x) != 0: return jn(n, x) / hn(n, x)
    else: return 0 

def f_sigma(x):
    lmbd = c / x #Лямбда
    k = 2 * pi / lmbd #Волновое число
    kr = k*r
    Sum = 0+0j
    for n in range(1, 40):
        Sum += (-1 ** n) * (n+ 1/2) * (bn(n, kr) - an(n, kr))
    return ((lmbd**2) / pi) * (np.abs(Sum)**2)

sigma = [f_sigma(x) for x in f]

plt.plot(f/1e9, sigma)
plt.xlabel('f, ГГц')
plt.ylabel('$\sigma, м^2$')
plt.grid()
plt.show()

