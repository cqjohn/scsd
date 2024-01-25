#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 16:58:17 2024

@author: genah
"""


#import libraries 
import matplotlib.pyplot as plt
import numpy as np 
from mpl_toolkits import mplot3d


def cosh_2d(x,y, xc, yc, T):
    return (np.cosh((x-xc)/T)**(-2))*(np.cosh((y-yc)/T)**(-2))

def add_source(h,w, xc, yc, T, data):  
    Y, X = np.ogrid[:h, :w]
    #dist_from_center = np.sqrt((X - mu[0])**2 + (Y-mu[1])**2)
    #mask = dist_from_center <= (1.5*sig)
    #rep = np.random.rand(*mask.shape)*np.max(sig*2) + np.ones((4611,2570))*(sig+mu)
    rep = cosh_2d(X, Y, xc, yc,T)
    data +=  rep
    return data 

x = np.linspace(-100, 100, 1000)
yr = np.cosh((x-10)/4)**(-2)
yy=np.cosh((x-10)/20)**(-2)
plt.plot(x,yr,'red')
plt.plot(x,yy,'yellow')
plt.show()


fig = plt.figure()
ax = plt.axes(projection='3d')



x = np.linspace(-2, 2, 1000)
y = np.linspace(-2, 2, 1000)
X, Y = np.meshgrid(x, y)
Z = cosh_2d(X, Y,0,0,2)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 50, cmap='binary')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');
plt.title("Cosh Broadening")
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
    Z = add_source(pd,pd, trslin[1][i],trslin[0][i],14, Z)  

contourplot = plt.contourf(X,Y,Z, cmap=plt.cm.bone,
      origin='lower')
ax.set_aspect('equal')
cbar = plt.colorbar(contourplot)
plt.show()


