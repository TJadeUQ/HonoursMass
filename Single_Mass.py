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
Mass_gal = np.array([df.MassM_0])
print(Mass_gal.shape)

#Convert redshift to comoving

#Constants
c = 3*10**5 #km/s
H_0 = 70 #km/s/Mpc
L_r0 = 0.75
L_m0 = 0.25
L_k0 = 0.044
L_l0 = 1

#Comoving
def Comoving(z_0,z_1):
    c = 3*10**5 #km/s
    H_0 = 70 #km/s/Mpc
    L_r0 = 0.75
    L_m0 = 0.25
    L_k0 = 0.044
    L_l0 = 1
    
    chi_z = lambda z: (c/H_0)*((L_k0 *(1+z)**2 + L_m0*(1+z)**3 + L_r0 * (1+z)**4 + L_l0))**(-1/2)

    Comov = quad(chi_z, z_0, z_1)
    return(Comov)

chi_1 = Comoving(0, min(z_red))
chi_2 = Comoving(0, max(z_red))
print('Comoving distances',chi_1,chi_2)

Co_diff = chi_2[0]-chi_1[0]
print('Co moving difference for volume',Co_diff,'Mpc')

#Calculate the density of the sample


#Volume in R_0
radius = Co_diff * np.tan(3*np.pi/180)
print('radius cone',radius)
V = (np.pi/3) * ((radius)**2) * Co_diff
print('Volume',V,'Mpc^3')

#Mass in M_0 (from spreadsheet)
M = sum(Mass_gal)
print('Mass total',M, 'M_0')

#density in M_0/R_0^3, maths checks out
rho = M/V
print(rho, 'M_0/Mpc^3')


#Want to calculate the numerical value of the velocity of each object
#Smoothing function
def SmthFunc(Ri, R, rs):
    
    Smth = np.array([])
    for j in range(len(Ri)):
        if Ri[j] != R:
            if la.norm(Ri[j]-R)<rs:
                Sm_Func = (abs(Ri[j]-R)**3)/rs**3   
            else:
                Sm_Func = 1
            Smth = np.append(Smth, Sm_Func)
        else:
            pass
    return(Smth)


def Three_r(Ri, R):
    r_cubed_term = np.array([])
    for y in Ri:
        Diff = (y - R)/((abs(y-R))**3)
        r_cubed_term = np.append(r_cubed_term, Diff)
    return(r_cubed_term)



def Velocity(Smth, r_cubed_term, W, R, density, beta_0):
    Trial_Sum = np.array([])
    for j in range(len(W)):
        Int = Smth[j] * W[j] * r_cubed_term[j]
        Trial_Sum = np.append(Trial_Sum, Int)
    
    Total = beta_0*(1/(4*3.14*density)*(np.nansum(Trial_Sum)) + R/3)
    
    return(Total)

#Calculate radial coordinates
"""
r_i_in = np.array([])
for i,j,k in zip(x_c,y_c,z_c):
    r = (i**2+j**2+k**2)**(1/2)
    r_i_in = np.append(r_i_in,r)
"""
    
r_i = np.array([])
error = np.array([])
for i in z_red:
    R = Comoving(0,i)
    r_i = np.append(r_i, R[0])
    error = np.append(error, R[1])
print('Comoving distance',r_i.shape,'Mpc/h') 


#Difference between radial and actual comoving distance
diff = np.array([])
for i, j in zip(r_i, d_c):
    d = i - j
    diff = np.append(diff,d)
print(diff)
#The difference between the actual and calculated values is a problem - the integration
#is probably the casue of this

#Run internal

Smooth = SmthFunc(r_i, r_i[0], min(r_i)/2)
print('Smooth',Smooth)
Pray = Three_r(r_i,r_i[0])
print('Cubed term',Pray.shape)
Vel = Velocity(Smooth, Pray, Mass_gal, r_i[0], rho, 0.5)
print('Velocity',Vel, 'MpcH_0/h', 'Assuming H_0 = 70km/s and h = 0.7, the velocity becomes', 10*Vel,'km/s')

#Compare the error
plt.errorbar([d_c],[r_i], yerr= [error], fmt = 'b', label = 'Comoving distance calclated using quad')
plt.errorbar([d_c],[d_c], yerr = 0, fmt = 'r', label = 'Comoving distance from MiceCat')
#plt.scatter([d_c],[r_i_in], label = 'Comoving distance from x,y,z coordinates')
plt.legend()
plt.xlabel('Comoving distance (Mpc/h)')
plt.ylabel('Comoving distance (Mpc/h)')

plt.show()


