import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
import scipy
from scipy.integrate import quad
import itertools as it
import pandas as pd
import cProfile
import re
from scipy import linalg as la
import array
from sympy import *
from mpl_toolkits.mplot3d import Axes3D
from collections import OrderedDict
from math import floor


df = pd.read_csv('Noice.csv')
z_red = df.z
x_c = df.x_c
y_c = df.y_c
z_c = df.z_c
vx = df.vx
vy = df.vy
vz = df.vz
vr = df.vr
log_m = df.log_m
d_c = df.d_c
ra = df.ra
Mass_gal = df.MassM_0

#Constants
c = 3*10**5 #km/s
H_0 = 70 #km/s/Mpc
L_r0 = 0.75
L_m0 = 0.25
L_k0 = 0.044
L_l0 = 1


#Create data array
data_array = np.array([[x_c], [y_c], [z_c]]).T

data_array_M = np.array([x_c,y_c,z_c, Mass_gal]).T
print(data_array_M)
data_array_1 = np.array([x_c, y_c, z_c]).T
#Create 3D histogram
hist, edges = np.histogramdd(data_array_1, bins = (3,3,3))
print(edges)
edges = np.array(edges)
print(edges)
#Use these edges to define data_array_M structure

#Create density field
def Mass_array(x_bin_1, x_bin_2, y_bin_1, y_bin_2, z_bin_1, z_bin_2, data):
    M = np.array([])
    for i in range(len(data)):
        if x_bin_1< data[i,0] <= x_bin_2 and y_bin_1 < data[i,1] <= y_bin_2 and z_bin_1 < data[i,2] <= z_bin_2:
            M = np.append(M, data[i,3])
            print(data[i,3])
        else:
            pass
    return(M)
print('Edge length',len(edges))
Mass = np.array([])
for K in range(len(edges)-1):
    for J in range(len(edges)-1):
        for I in range(len(edges)-1):
            print(I,J,K)
            F = Mass_array(edges[0,K], edges[0,K+1], edges[1,J], edges[1,J+1], edges[2, I], edges[2,I+1], data_array_M)
            print('Mass',F)
            Mass = np.append(Mass, sum(F))


print('FINAL',Mass)


#Size of bins
def Size(bins):
    Volume = np.array([])
    for I in range(len(bins)-1):
        for J in range(len(bins)-1):
            for K in range(len(bins)-1):
                V = ((bins[0,I+1]-bins[0,I])*(bins[1,J+1]-bins[1,J])*(bins[2,K+1]-bins[2,K]))
                Volume = np.append(Volume, V)
    return(Volume)
                
        
Cubes = Size(edges)
print(Cubes)

Density = np.array([])
#Density of bins
for Cube_Vol, Matter in zip(Cubes, Mass):
    Den = Matter/Cube_Vol
    Density = np.append(Density, Den)
print(Density)

#Overdenisity
#Total density
M_Tot = sum(Mass_gal)
V_Tot = sum(Cubes)
D_Tot = M_Tot/V_Tot

Over_Den = np.array([])
for den in Density:
    rho = (den - D_Tot) / D_Tot
    Over_Den = np.append(Over_Den, rho)
print(Over_Den)

Over_Den = np.reshape(Over_Den, np.shape(edges))
print(Over_Den)
# Smooth gridded distribution
def dosmooth(scale,datgrid,lx,ly,lz):
  print('\nSmoothing grid with scale =',scale,'...')
  nx,ny,nz = datgrid.shape[0],datgrid.shape[1],datgrid.shape[2]
  norm = np.sum(datgrid)
  x = lx*np.fft.fftfreq(nx)
  y = ly*np.fft.fftfreq(ny)
  z = lz*np.fft.fftfreq(nz)
  rsqgrid = x[:,np.newaxis,np.newaxis]**2 + y[np.newaxis,:,np.newaxis]**2 + z[np.newaxis,np.newaxis,:]**2
  kerngrid = np.exp(-rsqgrid/(2.*(scale**2)))
  datspec = np.fft.rfftn(datgrid)
  kernspec = np.fft.rfftn(kerngrid)
  datgrid = np.fft.irfftn(datspec*kernspec)
  datgrid *= norm/np.sum(datgrid)
 
  return datgrid
    
Smooth = dosmooth(1, Over_Den, len(edges[0]), len(edges[1]), len(edges[2]))
print(Smooth)

#Only outputting 8 masses - need to check the bins (works for (2,2,2) bins).
