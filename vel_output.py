import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import linalg as la
import itertools as it
import pandas as pd
import cProfile
import re

df = pd.read_csv('Noice.csv')
z = df.z
x_c = df.x_c
y_c = df.y_c
z_c = df.z_c
vx = df.vx
vy = df.vy
vz = df.vz
vr = df.vr
log_m = df.log_m
Mass_gal = df.MassM_0
print(Mass_gal.shape)
#define each of the masses
w_i = np.array([6.17*10**11, 1.99*10**12,2.43*10**11])
z_i = np.array([1.408,1.4158,1.4031])


#Convert redshift to comoving

#Constants
c = 3*10**8
H_0 = 0.7*100
L_r0 = 0.75
L_m0 = 0.25
L_k0 = 0.044
L_l0 = 1

#Comoving
def Comoving(Z):
    Dist = np.array([])
    for z in Z:
        D_c = c/H_0*(L_r0/5 * (1+z)**5 + L_m0/4 * (1+z)**4 + L_k0/3 *(1+z)**3 + L_l0 *z)
        Dist = np.append(Dist, D_c)
    return(Dist)
r_i = Comoving(z_i)
print(r_i)

#Want to calculate the numerical value of the velocity of each object
#Smoothing function
def SmthFunc(Ri, R, rs):
    Smth = np.array([])
    for i in range(len(R)):
        for j in range(len(Ri)):
            if la.norm(Ri[j]-R[i])<rs:
                Sm_Func = (la.norm(Ri[j]-R[i])**3)/rs**3
                
            else:
                Sm_Func = 1
            Smth = np.append(Smth, Sm_Func)
    Smth = np.reshape(Smth,(len(R),len(R)))
    return(Smth)


def Three_r(Ri, R):
    r_cubed_term = np.array([])
    for x in R:
        for y in Ri:
            Diff = (y - x)/(y-x)**3
            r_cubed_term = np.append(r_cubed_term, Diff)
    r_cubed_term = np.reshape(r_cubed_term, (len(Ri),len(Ri)))
    return(r_cubed_term)



def Velocity(Smth, r_cubed_term, W, R, rho, beta_0):
    Trial_Sum = np.array([])
    for j in range(len(R)):
        for i in range(len(R)):
            Int = Smth[j,i] * W[i] * r_cubed_term[j,i]
            Trial_Sum = np.append(Trial_Sum, Int)

    Trial_Sum = np.reshape(Trial_Sum, (len(R),len(R)))
    for j in range(len(R)):
        Total = np.array([])
        
        for i in range(len(R)):
            Tot = beta_0*(1/(4*3.14*rho)*(np.nansum(Trial_Sum)) + R[i]/3)
            Total = np.append(Total,('Mass', i+1, 'has velocity',  Tot))
    return(Total)


#Run internal

Smooth = SmthFunc(r_i, r_i, 500)
print(Smooth)
Pray = Three_r(r_i,r_i)
print(Pray)
Vel = Velocity(Smooth, Pray, w_i, r_i, 1,1)
print(Vel)



