# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 10:23:19 2014

@author: petar
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# Definiranje gaussovog impulsa koji propagira po vodu
def gauss(t):
    tu = 0.1
    sigma = 0.05
    gauss = np.exp(-(tu-t)**2/(2.*sigma**2))
    return gauss


print('FDTD simulacija putnog vala na idealnom vodu')

# Ulazni podaci:
L = 1.; T = 3.4
dx = 0.001; dt = 0.0005
v = 1.
fi1 = 0.; fi2 = 0.

# opcija = 1 <=> homogeni vod
# opcija = 2 <=> prijelaz zracni vod - kabel
# opcija = 3 <=> kabel umetnut u zracni vod
# opcija = 4 <=> prijelaz kabel u zracni vod
opcija = 2

print('Racunam ...')

Nx = int(L/dx)
Nt = int(T/dt)
r = np.zeros(Nx+1)
if opcija == 1:
    # Homogeni vod (bez promjene valne imped.)
    for i in range(Nx+1):    
        r[i] = ((v*dt)/dx)**2
elif opcija == 2:
    # Nehomogeni vod
    # Val iz sredstva vece valne impedancije prodire u sredstvo manje valne impedancije
    for i in range(Nx+1):
        if i <= Nx/2:
            r[i] = ((v*dt)/dx)**2
        else:
            r[i] = (((v/2.)*dt)/dx)**2
elif opcija == 3:
    # Nehomogeni vod
    # Dvije promjene sredstva (na svakoj trecini)
    # Srednji dio je manje valne impedancije od krajnjih dijelova
    # (npr. kabelski vod umetnut u dalekovod)
    for i in range(Nx+1):
        if i <= Nx/3:
            r[i] = ((v*dt)/dx)**2
        elif i >= 2*Nx/3:
            r[i] = ((v*dt)/dx)**2
        else:
            r[i] = (((v/2.)*dt)/dx)**2
elif opcija == 4:
    # Nehomogeni vod
    # Val iz sredstva manje valne impedancije prodire u sredstvo vece valne impedancije
    for i in range(Nx+1):
        if i <= Nx/2:
            r[i] = (((v/2.)*dt)/dx)**2
        else:
            r[i] = ((v*dt)/dx)**2

# Inicijalizacija, definiranje rubnih i pocetnih uvjeta
fi = np.zeros(Nx*Nt).reshape(Nx,Nt)

for j in range(1,Nt):
    fi[0,j] = fi1
    fi[Nx-1,j] = fi2

for j in range(Nt):
    ti = float(j)*dt
    fi[0,j] = gauss(ti)

for i in range(1,Nx-1):
    fi[i,1] = (1.-r[i])*fi[i,0] + (r[i]/2.)*(fi[i-1,0]+fi[i+1,0])

# FDTD simulacija
for j in range(1,Nt-1):
	for i in range(1,Nx-1):
		fi[i,j+1] = (2.*(1.-r[i]))*fi[i,j] + r[i]*(fi[i+1,j]+fi[i-1,j]) - fi[i,j-1]

print('Animacija ...')

# --------------------------------
# Animacija koristenjem matplotlib
# -------------------------------- 
def init():
    line.set_ydata(np.ma.array(x, mask=True))
    return line,
def animate(i):
    line.set_ydata(fi[:,i])  # update the data
    return line,

if opcija == 1:
    # Homogeni vod
    fig, ax = plt.subplots()
    x = np.linspace(0., L, Nx)  # x-axis
    line, = ax.plot(x, fi[:,0])
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-1., 1.)
    ax.grid(True)
    ani = animation.FuncAnimation(fig, animate, np.arange(0, Nt, 20, dtype=int), 
                                  init_func=init, blit=True)
    plt.show()

elif opcija == 2:
    # Nehomogeni vod (2 voda)
    fig, ax = plt.subplots()
    x0 = np.zeros(2)+0.5
    y0 = np.linspace(-1., 1., 2)
    x = np.linspace(0., L, Nx)  # x-axis
    line, = ax.plot(x, fi[:,0])
    ax.plot(x0, y0, 'k--')
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-1., 1.)
    ax.grid(True)
    ani = animation.FuncAnimation(fig, animate, np.arange(0, Nt, 20, dtype=int), 
                                  init_func=init, blit=True)
    plt.show()

elif opcija == 3:
    # Nehomogeni vod (3 voda)
    fig, ax = plt.subplots()
    x0 = np.zeros(2)+0.33
    y0 = np.linspace(-1., 1., 2)
    x1 = np.zeros(2)+0.66
    y1 = np.linspace(-1., 1., 2)
    x = np.linspace(0., L, Nx)  # x-axis
    line, = ax.plot(x, fi[:,0])
    ax.plot(x0, y0, 'k--')
    ax.plot(x1, y1, 'k--')
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-1., 1.)
    ax.grid(True)
    ani = animation.FuncAnimation(fig, animate, np.arange(0, Nt, 20, dtype=int), 
                                  init_func=init, blit=True)
    plt.show()

elif opcija == 4:
    # Nehomogeni vod (2 voda)
    fig, ax = plt.subplots()
    x0 = np.zeros(2)+0.5
    y0 = np.linspace(-1.5, 1.5, 2)
    x = np.linspace(0., L, Nx)  # x-axis
    line, = ax.plot(x, fi[:,0])
    ax.plot(x0, y0, 'k--')
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-1.5, 1.5)
    ax.grid(True)
    ani = animation.FuncAnimation(fig, animate, np.arange(0, Nt, 20, dtype=int), 
                                  init_func=init, blit=True)
    plt.show()
    
print('Kraj programa.')
