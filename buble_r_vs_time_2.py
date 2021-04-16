# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:41:59 2020

@author: Alex Lee
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline


A = 17.7
# B = 0.0036
B = 0.0086 ## LJ fluid

def drdt(w,t, ts, v1, rho1, gamma):
    r = w
    drdt = -(A**2*np.sqrt(t-ts)/B + 2*v1/r) + np.sqrt( A**2-2*gamma/(rho1*r) + (2*v1/r + A**2*np.sqrt( t-ts )/B)**2 )
    return drdt

def d2rdt2(r,drdt, t, ts, v1, rho1, gamma):

    term1 = -A**2/(2*B*np.sqrt(t-ts))
    term2 = 2*v1*drdt/r**2
    term3n = 2*gamma*drdt/(rho1*r**2) + 2*(A**2*np.sqrt(t-ts)/B+2*v1/r)*(-term1-term2)
    term3d = 2*np.sqrt(A**2+(A**2*np.sqrt(t-ts)/B+2*v1/r)**2 - 2*gamma/(rho1*r))
    d2rdt = term1 + term2 + term3n/term3d
    return d2rdt

mu_m = 1e-6
    
def find_ic(lcyl, Edep):
    ts = 0.082*lcyl
    rs = 0.38*lcyl + 217*Edep/lcyl

    return ts, rs

## 5.6mev 42600 nm
## 6.1mev 48300 nm
## 7.8mev 70000 nm
    
## 5.5kev 21.8 nm
## 50.5kev 187.7 nm
## 500.5kev 1740 nm
    
t0, r0 = find_ic(70000, 7800) ## nm & keV

t0 = t0*1e-9 ## in s
r0 = r0*1e-9 ## in m
# t = np.linspace(t0, 1, 1000000)
t = np.logspace(np.log10(t0), 0, 1000000)
# v1 = 0.141e-6
v1 = 18.3e-6    ##LJ fluid
rho1 = 1379
# gamma = 4.9*1e-3 
gamma = 3.667e-3 ##LJ fluid
ts = t[0]

sol = odeint(drdt, r0, t, args=(ts, v1, rho1, gamma))
sol = sol.flatten()

rt = CubicSpline(t, sol)

plt.figure(dpi = 200)
plt.plot(t, rt(t))
plt.xscale('log')
plt.yscale('log')
plt.grid()


c1 = 333 ## m/s

# r = sol

# np.savetxt("5.6Mev.csv", (t ,sol), delimiter=",", fmt=('%s, %f'))

# drdt1 = drdt(rt(t), t, ts, v1, rho1, gamma)
# drdt_cs = CubicSpline(t, drdt1)

# plt.figure()
# plt.plot( t,drdt1)
# plt.xscale('log')
# plt.yscale('log')
# d2rdt2 = d2rdt2(rt(t), drdt1, t, ts, v1, rho1, gamma)

# plt.figure()
# plt.plot(t, drdt_cs(t,1))
# plt.xscale('log')
# plt.xlim(np.log10(t0), 1e-4)

# plt.figure()
# volum_pp = rt(t)**2*d2rdt2 + 2*rt(t)*drdt_cs(t)**2
# plt.plot(t, volum_pp)
# # plt.ylim(-1e-5, 1e-4)
# plt.xscale('log')
# plt.yscale('log')
# P = 4*np.pi*rho1/(c1)*(volum_pp)**2

# # dVdt, dt = dfdt(V, t)
# # dVdt2, dt2 =  dfdt(dVdt, dt)

# # P = rho1*dVdt2**2/(4*np.pi*c1)
# plt.figure()
# plt.plot(t, P)
# plt.xscale('log')
# plt.yscale('log')
# plt.grid()

# r_save = np.transpose(np.array([t, rt(t)]))
# files = open('r_7.8mev.npy', 'wb')
# np.save(files, r_save)




# v_save = np.transpose(np.array([t, drdt1]))
# files = open('6.1mev.npy', 'wb')
# np.save(files, v_save)

# read = open('r_5.5kev.npy', 'rb')
# data = np.load(read)
# plt.plot(data[:, 0], data[:, 1])
# plt.xscale('log')
# plt.yscale('log')