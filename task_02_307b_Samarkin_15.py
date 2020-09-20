from scipy.constants import c, pi
from scipy.special import spherical_jn as jn, spherical_yn as yn
import numpy as np
import pathlib as pth
import matplotlib as plt

#Тестовые переменные
D = 10e-3
fmin = 0.01e9
fmax = 40e9

f = np.linspace(fmin, fmax, 400)
def hn(n, x):
    return jn(n, x) + 1j*yn(n,x)
print(f)
