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
n_bins = 4
#Create data array

data_array = np.array([[x_c], [y_c], [z_c]]).T
print(data_array)
data_array_M = np.array([x_c,y_c,z_c, Mass_gal]).T

data_array_1 = np.array([x_c, y_c, z_c]).T
#Create 3D histogram
hist, edges = np.histogramdd(data_array_1, bins = (n_bins,n_bins,n_bins))
print(hist)
print(edges)
edges = np.array(edges)
print(len(edges), edges.shape[1])
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
            print(edges[0,K], edges[0,K+1], edges[1,J], edges[1,J+1], edges[2, I], edges[2,I+1])
            
            F = Mass_array(edges[0,K], edges[0,K+1], edges[1,J], edges[1,J+1], edges[2, I], edges[2,I+1], data_array_M)
            print('Sum',F)
            Mass = np.append(Mass, sum(F))

print('FINAL',Mass)

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

Over_Den = np.reshape(Over_Den, [(len(edges[1])-1),(len(edges[1])-1),(len(edges[1])-1)])
print('OVER',Over_Den)


# Smooth gridded distribution
def dosmooth(scale,datgrid,lx,ly,lz):
  #print('\nSmoothing grid with scale =',scale,'...')
  nx,ny,nz = datgrid.shape[0],datgrid.shape[1],datgrid.shape[2]
  norm = np.sum(datgrid)
  x = lx*np.fft.fftfreq(nx)
  y = ly*np.fft.fftfreq(ny)
  z = lz*np.fft.fftfreq(nz)
  rsqgrid = x[:,np.newaxis,np.newaxis]**2 + y[np.newaxis,:,np.newaxis]**2 + z[np.newaxis,np.newaxis,:]**2
  kerngrid = np.exp(-rsqgrid/(2.*(scale**2))) #Gaussian smoothing kernel
  datspec = np.fft.rfftn(datgrid)
  kernspec = np.fft.rfftn(kerngrid) #Fourier transform of smoothed kernel
  multiply = datspec*kernspec
  datgrid = np.fft.irfftn(multiply)
  datgrid *= norm/np.sum(datgrid)
  return datgrid, x, y, z


#Smoothed 3D density field  
Smooth, k_space_x, k_space_y, k_space_z = dosmooth(1, Over_Den, n_bins, n_bins, n_bins)
print('K-Space',k_space_x)
print('Smooth', Smooth, Smooth.shape)


#Velocity field


def Vel_General(OverDen_Gen, k_space, k_mag, a, f, H):
    v_k = lambda k: - complex(0,1) * a * f * H * k / k_mag**2 * OverDen_Gen #Equation 5
    v_real = np.fft.ifftn(v_k(k_space)) #Fourier transform velocity into real space
    return v_real

#Break down into components 

#Calculate the magnitude of wavevector
k_mag_x = la.norm(k_space_x)
k_mag_y = la.norm(k_space_y)
k_mag_z = la.norm(k_space_z)

#Calculate each component       
vel_x = Vel_General(Smooth,k_space_x, k_mag_x, 1, 1, 70)
vel_y = Vel_General(Smooth, k_space_y, k_mag_y, 1, 1, 70)
vel_z = Vel_General(Smooth, k_space_z, k_mag_z, 1, 1, 70)
print('Vel',vel_x, vel_y, vel_z)

#Plot the actual and calculated velocities
fig = plt.figure()
ax = plt.axes(projection = '3d')

ax.scatter(vx, vy, vz, alpha = 0.3)
ax.set_xlabel('X-velocity (Mpc/s?)')
ax.set_ylabel('Y-velocity (Mpc/s?)')
ax.set_zlabel('Z-velocity (Mpc/s?)')

plt.show()

ax.scatter(vel_x, vel_y, vel_z, alpha = 0.3)
ax.set_xlabel('X-velocity (Mpc/s?)')
ax.set_ylabel('Y-velocity (Mpc/s?)')
ax.set_zlabel('Z-velocity (Mpc/s?)')
plt.show()

"""
#Mid point of bins
mid = np.array([])

for K in range(len(edges[1])-1):
    for J in range(len(edges[1])-1):
        for I in range(len(edges[1])-1):
            print(edges[0,K], edges[0,K+1], edges[1,J], edges[1,J+1], edges[2, I], edges[2,I+1])
            mx = (edges[0,I+1]+edges[0,I])/2
            my = (edges[1,J+1]+edges[1,J])/2
            mz = (edges[2,K+1]+edges[2,K])/2
            mid = np.append(mid, np.array([mx,my,mz]))
            print('Midpoint',mid)

mid = mid.reshape(n_bins**(3), 3)
print(mid, mid.shape)

print(mid[:,1])

scatter = ax.scatter(mid[:,0], mid[:,1], mid[:,2], c = Smooth_Colour, cmap = 'viridis_r', alpha = 0.3);
ax.set_xlabel('X-coordinates (Mpc)')
ax.set_ylabel('Y-coordinates (Mpc)')
ax.set_zlabel('Z-coordinates (Mpc)')
ax.set_xlim(0,max(x_c))
ax.set_ylim(0,max(y_c))
ax.set_zlim(0,max(z_c))

fig.colorbar(scatter)
plt.show()
"""
