"""""""""""""""""""""""""""""""""""""""""
Tyler Philp
University of Queensland
2018
Honours Project
"""""""""""""""""""""""""""""""""""""""""
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



#Want to split up df by the bins into arrays

#Create array?
data_array = np.array([z_red,Mass_gal])
"""
data_dict = {}
for i,j in zip(z_red, Mass_gal):
    data_dict[i] = [j]
print(data_dict)
"""
#print(data_array[:,0])


#Comoving distance between two redshifts
def Comoving(z_0,z_1):
    c = 3*10**5 #km/s
    H_0 = 70 #km/s/Mpc
    L_r0 = 0.75
    L_m0 = 0.25
    L_k0 = 0.044
    L_l0 = 1
    
    chi_z = lambda z: (c/H_0)*((L_k0 *(1+z)**2 + L_m0*(1+z)**3 + L_r0 * (1+z)**4 + L_l0))**(-1/2)

    Comov = quad(chi_z, z_0, z_1)
    #print(Comov)
    return(Comov)


#For subsequent Bins, will need volume of frustrum

#print('density', rho, 'M_0/mpc^3')

#Want to relate the redshifts - take the halves of eah violume (with z2)
def Vol(z_1,z_2):
    V = np.pi/3 *(z_2-z_1) *((z_2*np.tan(3*np.pi/180))**2 + z_1*np.tan(3*np.pi/180)*z_2*np.tan(3*np.pi/180) + (z_1*np.tan(3*np.pi/180))**2)
    return(V)

def Same_Vol(z_1, z_3):
    z_2 = ((z_1**3 + z_3**3)/2)**(1/3)
    #Volume for reshift frustrum
    Vol_1 = Vol(z_1,z_2)
    Vol_2 = Vol(z_2,z_3)
    
    if round(Vol_1) == round(Vol_2):
        print('True')
    else:
        print('False')
    return(z_1,z_2,z_3)

Vol_full = Same_Vol(min(data_array[0,:]),max(data_array[0,:]))

#Find some volumes
#Volume A
z_1_A = Vol_full[0]
z_3_A = Vol_full[1]

Vol_A = Same_Vol(z_1_A, z_3_A)

#Volume B
z_1_B = Vol_full[1]
z_3_B = Vol_full[2]

Vol_B = Same_Vol(z_1_B, z_3_B)
#########################
#Volume C
z_1_C = z_1_A
z_3_C = Vol_A[1]

Vol_C = Same_Vol(z_1_C, z_3_C)

#Volume D
z_1_D = Vol_A[1]
z_3_D = Vol_A[2]

Vol_D = Same_Vol(z_1_D, z_3_D)

#Volume E
z_1_E = Vol_A[2]
z_3_E = Vol_B[1]

Vol_E = Same_Vol(z_1_E, z_3_E)

#Volume F
z_1_F = Vol_B[1]
z_3_F = Vol_B[2]

Vol_F = Same_Vol(z_1_F, z_3_F)

#Create Bins
Bins = np.array([Vol_A, Vol_B, Vol_C, Vol_D, Vol_E, Vol_F])

Bins = np.reshape(Bins, (1,-1))

Bins = np.unique(Bins)

print('RED BINS',Bins)

#Plot with new bins
red_bins, bin_edges = np.histogram(z_red, bins = Bins)

plt.hist(z_red, bins = Bins)
plt.xlabel('Redshift bins for constant volume')
plt.ylabel('Galaxy Count')
plt.show()

#########################################
#Density of each bin
#########################################

#Need to find the mass in each bin
data_dict = {}
for i,z in zip(Mass_gal,z_red):
    data_dict[z] = [i]        


def Mass_array(z_1, z_2, dictionary):
    d = np.array([])
    
    for key in dictionary:
        if z_1 < key <= z_2:
            L = dictionary[key]
            d = np.append(d, L)
    return(d)


d1 = sum(Mass_array(Bins[0], Bins[1], data_dict))
d2 = sum(Mass_array(Bins[1], Bins[2], data_dict))
d3 = sum(Mass_array(Bins[2], Bins[3], data_dict))
d4 = sum(Mass_array(Bins[3], Bins[4], data_dict))
d5 = sum(Mass_array(Bins[4], Bins[5], data_dict))
d6 = sum(Mass_array(Bins[5], Bins[6], data_dict))
d7 = sum(Mass_array(Bins[6], Bins[7], data_dict))
d8 = sum(Mass_array(Bins[7], Bins[8], data_dict))

M = np.array([d1,d2,d3,d4,d5,d6,d7,d8])

#Density

Volume = np.array([])

for i in range(len(Bins)-1):
    V = Vol(Bins[i], Bins[i+1])
    Volume = np.append(Volume, V)
    

rho = np.array([])
for L, A in zip(Volume, M):
    den = A/L
    rho = np.append(rho, den)
print('RHO', rho)

#Overdensity
Total_Mass = sum(Mass_array(Bins[0], Bins[8], data_dict))
Total_Vol = Vol(Bins[0], Bins[8])
Total_Den = Total_Mass/Total_Vol
print('Tota',Total_Den)

Over_Den = np.array([])
for i in range(len(rho)):
    O = rho[i] - Total_Den
    Over_Den = np.append(Over_Den, O)
print('OVER',Over_Den)

#########################################
#FFT!
#########################################

#Gaussian
def Gaussian(x,sigma, mu):
    Gauss = 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-1/2 * ((x-mu)/sigma)**2)
    return(Gauss)

def variance(X,mu,N):
    diff = np.array([])
    for i in X:
        Summa = (i-mu)
        diff = np.append(diff, Summa)
    diff2 = np.array([])
    for j in diff:
        diff2 = np.append(diff2, j**2)
    Sigma2 = sum(diff2)/N
    
    return(Sigma2)

#Breakup data into bins and multiply by gaussian
x = data_array[0,:]

s= (variance(x, np.mean(x), len(x)))**(1/2)
print('Variance', s)


F = np.array([])
G = np.array([])
for i in x:
    if Bins[0]<i<Bins[1]:
        value = Over_Den[0]
        F = np.append(F, value)
        Ga = Gaussian(i, s, np.mean([bin_edges[0], bin_edges[1]]))
        G = np.append(G, Ga)
    elif Bins[1]<i<Bins[2]:
        value = Over_Den[1]
        F = np.append(F, value)
        Ga = Gaussian(i, s, np.mean([bin_edges[1], bin_edges[2]]))
        G = np.append(G, Ga)
    elif Bins[2]<i<Bins[3]:
        value = Over_Den[2]
        F = np.append(F, value)
        Ga = Gaussian(i, s, np.mean([bin_edges[2], bin_edges[3]]))
        G = np.append(G, Ga)
    elif Bins[3]<i<Bins[4]:
        value = Over_Den[3]
        F = np.append(F, value)
        Ga = Gaussian(i, s, np.mean([bin_edges[3], bin_edges[4]]))
        G = np.append(G, Ga)
    elif Bins[4]<i<Bins[5]:
        value = Over_Den[4]
        F = np.append(F, value)
        Ga = Gaussian(i, s, np.mean([bin_edges[4], bin_edges[5]]))
        G = np.append(G, Ga)
    elif Bins[5]<i<Bins[6]:
        value = Over_Den[5]
        F = np.append(F, value)
        Ga = Gaussian(i, s, np.mean([bin_edges[5], bin_edges[6]]))
        G = np.append(G, Ga)
    elif Bins[6]<i<Bins[7]:
        value = Over_Den[6]
        F = np.append(F, value)
        Ga = Gaussian(i, s, np.mean([bin_edges[6], bin_edges[7]]))
        G = np.append(G, Ga)
    elif Bins[7]<i<Bins[8]:
        value = Over_Den[7]
        F = np.append(F, value)
        Ga = Gaussian(i, s, np.mean([bin_edges[7], bin_edges[8]]))
        G = np.append(G, Ga)

    else:
        value = 0
        F = np.append(F, value)
        G = np.append(G, value)



H = np.array([])
for i,j in zip(G,F):
    h = i*j
    H = np.append(H,h)


#Fourier transform of Gaussian and overdensity - produce delta_m
H_fft = np.fft.fft(H)


#Fourier and real
#plt.plot(x,H_fft)
plt.plot(x,F)
plt.plot(x,H)
#plt.plot(x,H_ifft)
plt.xlabel('Redshift')
#plt.ylabel('Fourier transform and actual count')
plt.gca().legend(('Overdensity','Fourier transform'))
plt.show()

#Calculate the velocity in k space
#v(k) = -iaHfk/k^2 theta(k)
#delta_g = b * delta_m

#Let's just do the radial parts
#Comoving distance?
r = np.array([])
for i in range(len(z_red)):
    C = Comoving(0, z_red[i])
    r = np.append(r,C[0])

k_r = np.fft.fft(r)
print(k_r)
#mod_k_r = np.linalg.norm(k_r)

a = 1
f = 1

#MAIN EVENT: CALCULATE THE VELOCITY IN K SPACE
v_r_k = np.array([])
for K in range(len(r)):
    v = -(a*H_0*f)*(k_r[K]/(k_r[K])**2) * H_fft[K]
    v_Im = np.imag(v)
    v_r_k = np.append(v_r_k, v_Im)
print('Radial Fourier vel',v_r_k)
#Transform velocity to real space
v_r_R = np.real(np.fft.fft(v_r_k))
print(v_r_R)



"""
k_x = np.fft.fft(x_c)
k_y = np.fft.fft(y_c)
k_z = np.fft.fft(z_c)

k = np.array([k_x,k_y,k_z])
k = np.reshape(k, (2721,3))

mod_k = np.array([])
for vec in k:
    modulus = np.linalg.norm(vec)
    mod_k = np.append(mod_k, modulus)
    

#Component velocities
a = 1
f = 1
v_x_k = np.array([])
for i in range(len(k_x)):
    v = -((a*H_0*f)*(k_x[i]/(mod_k[i])**2) * H_fft[i])
    v_x_k = np.append(v_x_k, v)

v_y_k = np.array([])
for j in range(len(k_y)):
    v = -((a*H_0*f)*(k_y[j]/(mod_k[j])**2) * H_fft[j])
    v_y_k = np.append(v_y_k, v)

v_z_k = np.array([])
for m in range(len(k_z)):
    v = -((a*H_0*f)*(k_z[m]/(mod_k[m])**2) * H_fft[m])
    v_z_k = np.append(v_z_k, v)

v_k = np.array([v_x_k, v_y_k, v_z_k])


#Transform back to real space
v_x_R = np.real(np.fft.fft(v_x_k)).tolist()

v_y_R = np.real(np.fft.fft(v_y_k)).tolist()

v_z_R = np.real(np.fft.fft(v_z_k)).tolist()

v_R = (v_x_R, v_y_R, v_z_R)

v_hist_R, bin_edges_R = np.histogramdd(v_R)

fig = plt.figure()
ax1 = fig.add_subplot(111, projection = '3d')

ax1.set_xlabel('X Velocity')
ax1.set_ylabel('Y Velocity')
ax1.set_zlabel('Z Velocity')

ax1.plot(v_x_R, v_y_R, v_z_R, 'k.', alpha = 0.1)
ax1.plot(vx,vy,vz, 'r.', alpha = 0.1)
"""
'''
print('HOLLA',bin_edges_R[0][:-1])
X, Y = np.meshgrid(bin_edges_R[0][:-1], bin_edges_R[1][:-1])

for ct in [0,1,2,3,4,5,6,7,8]:
    cs = ax1.contourf(X,Y,v_hist_R[ct], zdir = 'z', offset = bin_edges_R[2][ct], cmap = plt.cm.RdYlBu_r, alpha = 0.5)
plt.colorbar(cs)

ax1.set_xlim(min(v_x_R), max(v_x_R))
ax1.set_ylim(min(v_y_R), max(v_y_R))
ax1.set_zlim(min(v_z_R), max(v_z_R))
'''
plt.scatter(z_red, vr, alpha = 0.3)
plt.scatter(z_red, v_r_R, alpha = 0.3)

plt.show()


#Difference
difference = np.array([])
for R,D in zip(vr, v_r_R):
    d = R-D
    difference = np.append(difference, d)

plt.plot(z_red, difference)
plt.xlabel('Redshift')
plt.ylabel('Difference in peculiar velocity (real-Calculated)')
plt.show()
