# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 16:51:50 2020

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
    
def cal_vsrc(lcyl, Edep, filename): 
    v1 = 18.3e-6    ##LJ fluid
    rho1 = 1379
    # gamma = 4.9*1e-3 
    gamma = 3.667e-3 ##LJ fluid
    A = 17.7
    # B = 0.0036
    B = 0.0086 ## LJ fluid
    
    t0, r0 = find_ic(lcyl, Edep) ## nm & keV 
    
    t0 = t0*1e-9 ## in s
    r0 = r0*1e-9 ## in m 
    t =np.linspace(t0, 1e-2, 1000000)
    sol = odeint(drdt, r0, t, args=(t0, v1, rho1, gamma))
    sol = sol.flatten() 
    # print(sol)
    rt = CubicSpline(t, sol)
    c1 = 333 ## m/s
    drdt1 = drdt(rt(t), t, t0, v1, rho1, gamma)
    
    v_save = np.transpose(np.array([t, drdt1]))
    files = open(filename, 'wb')
    np.save(files, v_save)
    print(filename + ' saved')
    
    
cal_vsrc(42600, 5.6e3, '5.6mev.npy')    
cal_vsrc(48300, 6.1e3, '6.1mev.npy')    
cal_vsrc(70000, 7.8e3, '7.8mev.npy')    

cal_vsrc(21.8, 5.5, '5.5kev.npy')    
cal_vsrc( 187.7, 50.5,'50.5kev.npy')    
cal_vsrc( 1740, 500.5, '500.5kev.npy')    
    
    
def d2rdt2(r,drdt, t, ts, v1, rho1, gamma):

    term1 = -A**2/(2*B*np.sqrt(t-ts))
    term2 = 2*v1*drdt/r**2
    term3n = 2*gamma*drdt/(rho1*r**2) + 2*(A**2*np.sqrt(t-ts)/B+2*v1/r)*(-term1-term2)
    term3d = 2*np.sqrt(A**2+(A**2*np.sqrt(t-ts)/B+2*v1/r)**2 - 2*gamma/(rho1*r))
    d2rdt = term1 + term2 + term3n/term3d
    return d2rdt    
    
    
    
    
    