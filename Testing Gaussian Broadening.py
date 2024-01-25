#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 16:28:21 2024

@author: genah
"""

#import libraries 
import matplotlib.pyplot as plt
import numpy as np 



def gaussian(x, mu, sig, A):
    return A*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def gaussian_2d(x, y, mu, sig, A):
    mx, my = mu[0], mu[1] 
    sx, sy = sig, sig 
    #/ (2. * np.pi * sx * sy)
    return  A  * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.)))
      
    
def add_source(h,w, mu, sig, A, data):  
    Y, X = np.ogrid[:h, :w]
    #dist_from_center = np.sqrt((X - mu[0])**2 + (Y-mu[1])**2)
    #mask = dist_from_center <= (1.5*sig)
    #rep = np.random.rand(*mask.shape)*np.max(sig*2) + np.ones((4611,2570))*(sig+mu)
    rep = gaussian_2d(X, Y, mu, sig, A)
    data +=  rep
    return data 

# Data for a three-dimensional line
x = np.linspace(-50, 50, 1000)
y = np.linspace(-50, 50, 1000)
X, Y = np.meshgrid(x, y)
Z = gaussian_2d(X, Y, [10,10], 5, 10)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 50, cmap='binary')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');
plt.title("Gaussian Broadening")
plt.show()



pd = 1000
Vg1 = np.linspace(0,5,pd)#*1e-6
Vg2 = np.linspace(0,5,pd)#*1e-6 

X, Y = np.meshgrid(Vg1, Vg2)
fig, ax = plt.subplots()


Z = np.zeros([pd,pd])
#Z= np.random.rand(*Z.shape)*np.max(180) + np.ones((pd,pd))*100 #add white noise
Z[500:501] = np.ones(pd)


trslin = np.where(Z ==1)

#[[trslin[0][i],trslin[1][i]] for _ in range(0,len(trslin))]
for i in range(0,len(trslin[0])):  
    print(i,'/',len(trslin[0]))
    Z = add_source(pd,pd, [trslin[1][i],trslin[0][i]] , 20, 5, Z)  

contourplot = plt.contourf(X,Y,Z, cmap=plt.cm.bone,
      origin='lower')
ax.set_aspect('equal')
cbar = plt.colorbar(contourplot)
plt.show()
