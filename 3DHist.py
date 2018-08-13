import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
from scipy import linalg as la
import array
from mpl_toolkits import mplot3d
import pandas as pd


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
a = 1
f = 1
n_bins = 3
#Create data array

data_array = np.array([[x_c], [y_c], [z_c]]).T

data_array_M = np.array([x_c,y_c,z_c, Mass_gal]).T

data_array_1 = np.array([x_c, y_c, z_c]).T
#Create 3D histogram
hist, edges = np.histogramdd(data_array_1, bins = (n_bins,n_bins,n_bins))
edges = np.array(edges)
print(edges)


#Use these edges to define data_array_M structure

#Create density field
def Mass_array(x_bin_1, x_bin_2, y_bin_1, y_bin_2, z_bin_1, z_bin_2, data):
    M = np.array([])
    for i in range(len(data)):
        if x_bin_1< data[i,0] <= x_bin_2 and y_bin_1 < data[i,1] <= y_bin_2 and z_bin_1 < data[i,2] <= z_bin_2:
            M = np.append(M, data[i,3])
        else:
            pass
    return(M)

Mass = np.array([])
for K in range(len(edges[1])-1):
    for J in range(len(edges[1])-1):
        for I in range(len(edges[1])-1):           
            F = Mass_array(edges[0,K], edges[0,K+1], edges[1,J], edges[1,J+1], edges[2, I], edges[2,I+1], data_array_M)
            Mass = np.append(Mass, sum(F))


#Mass = Mass_gal
#Size of bins
def Size(bins):
    Volume = np.array([])
    for I in range(len(bins[1])-1):
        for J in range(len(bins[1])-1):
            for K in range(len(bins[1])-1):
                V = ((bins[0,I+1]-bins[0,I])*(bins[1,J+1]-bins[1,J])*(bins[2,K+1]-bins[2,K]))
                Volume = np.append(Volume, V)
    return(Volume)
                
        
Cubes = Size(edges)

Density = np.array([])
#Density of bins
for Cube_Vol, Matter in zip(Cubes, Mass):
    Den = Matter/Cube_Vol
    Density = np.append(Density, Den)


#Overdenisity
#Total density
M_Tot = sum(Mass_gal)
V_Tot = sum(Cubes)
D_Tot = M_Tot/V_Tot

Over_Den = np.array([])
for den in Density:
    rho = (den - D_Tot) / D_Tot
    Over_Den = np.append(Over_Den, rho)


Over_Den = np.reshape(Over_Den, [(len(edges[1])-1),(len(edges[1])-1),(len(edges[1])-1)])



# Smooth gridded distribution - Chris Blake code
def dosmooth(scale,datgrid,lx,ly,lz,b):

  nx,ny,nz = datgrid.shape[0],datgrid.shape[1],datgrid.shape[2]
  norm = np.sum(datgrid)

  x = lx*np.fft.fftfreq(nx)
  y = ly*np.fft.fftfreq(ny)
  z = lz*np.fft.fftfreq(nz)


  rsqgrid = x[:,np.newaxis,np.newaxis]**2 + y[np.newaxis,:,np.newaxis]**2 + z[np.newaxis,np.newaxis,:]**2 #Transform into radial coordinates
  kerngrid = np.exp(-rsqgrid/(2.*(scale**2))) #Gaussian smoothing kernel
  
  datspec = np.fft.fftn(datgrid)
  kernspec = np.fft.fftn(kerngrid) #Fourier transform of smoothed kernel
  
  mass_k = datspec*kernspec
  mass_k = b * mass_k

  datgrid = np.fft.ifftn(mass_k)
  datgrid *= norm/np.sum(datgrid)
  print('DAT',datgrid, datgrid.shape)
  
  kx = 2.*np.pi*np.fft.fftfreq(nx,d=lx/nx)
  ky = 2.*np.pi*np.fft.fftfreq(ny,d=ly/ny)
  kz = 2.*np.pi*np.fft.fftfreq(nz,d=lz/nz)#[:nz/2+1]


  ktot = np.sqrt(kx[:,np.newaxis,np.newaxis]**2 + ky[np.newaxis,:,np.newaxis]**2 + kz[np.newaxis,np.newaxis,:]**2)

  v_ktot = - H_0 * ktot/ (la.norm(ktot)**2) * mass_k

  v_rtot = np.fft.ifft(v_ktot)
  #print(v_rtot, v_rtot.shape)
  return v_rtot, datgrid

velreal, smooth_r = dosmooth(1, Over_Den, edges[0,1]-edges[0,0], edges[1,1]-edges[1,0], edges[2,1]-edges[2,0], 1)


#Radial theory
Real_Vel = velreal.real
Im_Vel = velreal.imag


Amp = np.sqrt(Real_Vel**2 + Im_Vel**2)
theta = np.arctan(Im_Vel/Real_Vel)

x_coord = Amp * np.sin(theta) * np.cos(theta)
y_coord = Amp * np.sin(theta) * np.sin(theta)
z_coord = Amp * np.cos(theta)

coord_vel = np.array([x_coord, y_coord, z_coord])

#Mid point of bins
midx = np.array([])
midy = np.array([])
midz = np.array([])
for I in range(n_bins):
    mx = (edges[0,I+1]+edges[0,I])/2
    midx = np.append(midx, mx)
for J in range(n_bins):
    my = (edges[1,J+1]+edges[1,J])/2
    midy = np.append(midy, my)
for K in range(n_bins):
    mz = (edges[2,K+1]+edges[2,K])/2
    midz = np.append(midz, mz)
            

#Plot the density field and velocity field
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

ax.quiver(midx, midy, midz, x_coord, y_coord, z_coord, length = 1)
plt.show()
"""
X,Y = np.meshgrid(edges[0][:-1],edges[1][:-1])
#Still casting the imaginary component of the density field?
for ct in [0,1,2]:
    cs = ax.contourf(X,Y,smooth_r[ct], zdir = 'z', offset = edges[2][ct],
                      cmap = plt.cm.RdYlBu_r, alpha = 0.5)
plt.colorbar(cs)
ax.set_xlim(0,max(midx))
ax.set_ylim(0,max(midy))
ax.set_zlim(0,max(midz))
plt.show()
"""
