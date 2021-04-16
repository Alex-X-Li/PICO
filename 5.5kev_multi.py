# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 17:59:16 2020

@author: Alex Lee
"""

import numpy as np
import multiprocessing as mp
from scipy.interpolate import CubicSpline
# import matplotlib.pyplot as plt


energy = '5.5kev'

def simu(ran):
    off = ran*5e-3  ###################### in m the sourse off from center
    # Discretization
    c1=10   # Number of grid points per dominant wavelength
    c2=0.2  # CFL-Number
    ## r = 150 mm
    ## z = 1000 mm
    fmax = 300e3 ## 300kHz
    vc3f8 = 333 ##m/s
    wl = 333/fmax
    gp = wl/c1
    
    nx = int(200e-3/gp)
    ny = int(200e-3/gp)
    
    
    
    # nx=300 # Number of grid points in X
    # ny=300 # Number of grid points in Y
    T=2e-4     # Total propagation time
    
    # Source Signal
    f0 = fmax      # Center frequency Ricker-wavelet
    q0= 100       # Maximum amplitude Ricker-Wavelet
    
    
    # Receiver
    # xrec1=int(0.5*nx); yrec1=int(0.5*ny*1.25);  # Position Reciever 1 (in grid points)
    # xrec2=int(0.5*nx); yrec2=int(0.5*ny*1.5);  # Position Reciever 2 (in grid points)
    # xrec3=int(0.5*nx); yrec3=int(0.5*ny*1.75);# Position Reciever 3 (in grid points)
    
    # Velocity and density
    vc3f8 = 333
    dx = vc3f8/(fmax*c1)             # Spatial discretization (in m)
    dy = dx                         # Spatial discretization (in m)
    vglass = 5760  ## m/s  sqrt(73.1e9/2203)
    vpiezo = 3183 ## m/s  sqrt(76e9/7500)
    # vvac = 346 #m/s
    voil = 1437 ##   sqrt(1.8e9/870.93) ISO 32 mineral oil
    
    rho_gla = 2203  ## kg/m3
    rho_piezo = 7500  ## kg/m3
    rho_c3f8 = 1397
    # rho_vac = 0.005842 
    rho_oil = 870.93 ####  kg/m³  Multitherm PG-1 Oil
    
    glass_dis = 0.01/dx ## 1cm
    glass_thick = 0.005/dx ## 5mm
    piezo_dia = 0.003175/dy ## mm
    pizo_hight = 0.0084/dx ## mm
    
    modell_v = np.ones((ny,nx))*vc3f8
    rho = np.ones((ny,nx))*rho_c3f8
    
    xscr = int(0.5*nx )  # Source position (in grid points) in X
    yscr = int(0.5*ny + off/dx)  # Source position (in grid points) in Y
    
    
    ## glass
    modell_v[int(nx//2+glass_dis):int(nx//2+glass_dis+glass_thick), : ] = vglass
    ## vacuum
    modell_v[int(nx//2+glass_dis+glass_thick): nx, : ] = voil
    ## piezo
    modell_v[int(nx//2+glass_dis+glass_thick): int(nx//2+glass_dis+glass_thick+pizo_hight), int(ny//2-piezo_dia) : int(ny//2+piezo_dia)] = vpiezo
    
    rho[int(nx//2+glass_dis):int(nx//2+glass_dis+glass_thick), : ] = rho_gla
    ## vacuum
    rho[int(nx//2+glass_dis+glass_thick): nx, : ] = rho_oil
    ## piezo
    rho[int(nx//2+glass_dis+glass_thick): int(nx//2+glass_dis+glass_thick+pizo_hight), int(ny//2-piezo_dia) : int(ny//2+piezo_dia)] = rho_piezo
    ## matrix mean : [y, x]
    modell_v = np.transpose(modell_v)
    rho = np.transpose(rho)
    # plt.imshow(modell_v)
    
    vx=np.zeros(shape = (ny,nx))
    vy=np.zeros(shape = (ny,nx))
    p=np.zeros(shape = (ny,nx))
    vx_x=np.zeros(shape = (ny,nx))
    vy_y=np.zeros(shape = (ny,nx))
    p_x=np.zeros(shape = (ny,nx))
    p_y=np.zeros(shape = (ny,nx))
    
    # Calculate first Lame-Paramter
    l=rho * modell_v * modell_v
    
    cmin=min(modell_v.flatten())  # Lowest P-wave velocity
    cmax=max(modell_v.flatten())  # Highest P-wave velocity
    #fmax=2*f0                     # Maximum frequency
    dx=cmin/(fmax*c1)             # Spatial discretization (in m)
    dy=dx                         # Spatial discretization (in m)
    dt=dx/(cmax)*c2               # Temporal discretization (in s)  ### change from cmax to cin
    lampda_min=cmin/fmax          # Smallest wavelength
    
    # Output model parameter:
    print("Model size: x:",dx*nx,"in m, y:",dy*ny,"in m")
    print("Temporal discretization: ",dt," s")
    print("Spatial discretization: ",dx," m")
    print("Number of gridpoints per minimum wavelength: ",lampda_min/dx)
    
    
    x=np.arange(0,dx*nx,dx) # Space vector in X
    y=np.arange(0,dy*ny,dy) # Space vector in Y
    t=np.arange(0,T,dt)     # Time vector
    nt=np.size(t)           # Number of time steps
    
    print('Number of time step: ', nt)
    print('Grid size:' ,nx, "x", ny)
    
    
    filename = energy + '.npy'
    files = open(filename, 'rb')
    data = np.load(files)
    t = data[:, 0]
    vdata = data[:, 1]
    t = t-t[0]
    vscr = CubicSpline(t, vdata)
    
    
    
    # Init Seismograms
    # Seismogramm=np.zeros((3,nt)); # Three seismograms
    
     # Calculation of some coefficients
    i_dx=1.0/(dx)
    i_dy=1.0/(dy)
    c1=9.0/(8.0*dx)
    c2=1.0/(24.0*dx)
    c3=9.0/(8.0*dy)
    c4=1.0/(24.0*dy)
    c5=1.0/np.power(dx,3)   ### dx**(-3)
    c6=1.0/np.power(dy,3)   ### dy**(-3)
    c7=1.0/np.power(dx,2)
    c8=1.0/np.power(dy,2)
    c9=np.power(dt,3)/24.0
     # Prepare slicing parameter:
    kxM2=slice(5-2,nx-4-2)
    kxM1=slice(5-1,nx-4-1)
    kx=slice(5,nx-4)
    kxP1=slice(5+1,nx-4+1)
    kxP2=slice(5+2,nx-4+2)
    
    kyM2=slice(5-2,ny-4-2)
    kyM1=slice(5-1,ny-4-1)
    ky=slice(5,ny-4)
    kyP1=slice(5+1,ny-4+1)
    kyP2=slice(5+2,ny-4+2)
    
     # f_wr = open( energy + 'p_wave.npy', 'wb')
    
    savelen = int(nx//2+glass_dis+glass_thick+pizo_hight) - int(nx//2+glass_dis+glass_thick)
    saveline = np.zeros((savelen, nt)) ## 6000?
     # f_read = open( energy + '_'+ str(off) + '_p_wave.npy', 'rb')
    
     ## Time stepping
    print("Starting time stepping...")
     
    
    
    
        
        
        
        
    for n in range(2,nt):
            print(yscr, xscr)
             # Inject source wavelet
            vx[yscr,xscr+1] = vx[yscr,xscr+1] + vscr(n*dt)
            vx[yscr,xscr-1] = vx[yscr,xscr-1] - vscr(n*dt)
            vy[yscr+1,xscr] = vy[yscr+1,xscr] + vscr(n*dt)
            vy[yscr-1,xscr] = vy[yscr-1,xscr] - vscr(n*dt)
            
             # Update velocity
            p_x[ky,kx]=c1*(p[ky,kxP1]-p[ky,kx])-c2*(p[ky,kxP2]-p[ky,kxM1])
            p_y[ky,kx]=c3*(p[kyP1,kx]-p[ky,kx])-c4*(p[kyP2,kx]-p[kyM1,kx])
            
            vx=vx-dt/rho*p_x
            vy=vy-dt/rho*p_y
            
             # Update pressure
            vx_x[ky,kx]=c1*(vx[ky,kx]-vx[ky,kxM1])-c2*(vx[ky,kxP1]-vx[ky,kxM2])
            vy_y[ky,kx]=c3*(vy[ky,kx]-vy[kyM1,kx])-c4*(vy[kyP1,kx]-vy[kyM2,kx])
            
            p=p-l*dt*(vx_x+vy_y)
            
             # Save seismograms
             # Seismogramm[0,n]=p[yrec1,xrec1]
             # Seismogramm[1,n]=p[yrec2,xrec2]
             # Seismogramm[2,n]=p[yrec3,xrec3]
             # p_save[0].append(p)
             # if n-2 %  
             # np.save(f_wr, p)
            saveline[:, n] = p[ny//2, int(nx//2+glass_dis+glass_thick): int(nx//2+glass_dis+glass_thick+pizo_hight) ]
            
          
    print("Finished time stepping!")
     # np.save(energy + str(off) +'_last_step.npy', p)
     # print("Saved ", energy +'_last_step.npy' )
    
    np.save(energy +'_' + str(off) + "_pizoline.npy", saveline)
    print("Finished saving ", energy + str(off)  + "_pizoline.npy")
    



line = []

for i in range(0,2):
    off = 3 * np.random.ranf([1])
    proc = mp.Process(target=simu,args=(off,))
    line.append(proc)
    proc.start()
    
for i in line:
    i.join()



