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
vx_orig = df.vx
vy_orig = df.vy
vz_orig = df.vz
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
n_bins = 9
#Create data array
size = int(1000)
print(min(x_c),max(x_c),min(y_c),max(y_c),min(z_c),max(z_c))

#Extend the array to a predetermined field
""" 
x_c = np.pad(np.array(x_c), (size,size), 'linear_ramp', end_values = (-150,300))
y_c = np.pad(np.array(y_c), (size,size), 'linear_ramp', end_values = (-1000,5000))
z_c = np.pad(np.array(z_c), (size,size), 'linear_ramp', end_values = (-100,350))
Mass_gal = np.pad(np.array(Mass_gal), (size,size),'constant', constant_values = (0, 0))
"""  
Mass1= np.zeros(len(Mass_gal))
for i in range(len(x_c)):
    if 20.0<= x_c[i] <=130.0 and 500.0<=y_c[i]<=1000.0 and 20.0<= z_c[i] <=130.0:
        Mass1[i] = Mass_gal[i]
    else:
        Mass1[i] = 0

print(Mass1)
            
data_array = np.array([[x_c], [y_c], [z_c]]).T

data_array_M = np.array([x_c,y_c,z_c, Mass1]).T

data_array_1 = np.array([x_c, y_c, z_c]).T


#Create 3D histogram
hist, edges = np.histogramdd(data_array_1, bins = (n_bins,n_bins,n_bins))
edges = np.array(edges)
#print(edges)


#Use these edges to define data_array_M structure

#Create density field
def Mass_array(x_bin_1, x_bin_2, y_bin_1, y_bin_2, z_bin_1, z_bin_2, data):
    M = np.zeros(len(data))
    for i in range(len(data)):
        if x_bin_1< data[i,0] <= x_bin_2 and y_bin_1 < data[i,1] <= y_bin_2 and z_bin_1 < data[i,2] <= z_bin_2:
            M[i] = data[i,3]
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
    rho = den / D_Tot
    Over_Den = np.append(Over_Den, rho)


Over_Den = np.reshape(Over_Den, [(len(edges[1])-1),(len(edges[1])-1),(len(edges[1])-1)])

print("D_TOT", D_Tot)

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

  datgrid = np.fft.ifftn(datspec*kernspec)
  datgrid *= norm/np.sum(datgrid)

  mass_k = np.fft.fftn(datgrid)
 
  mass_k = b * mass_k

  kx = 2.*np.pi*np.fft.fftfreq(nx,d=lx/nx)
  ky = 2.*np.pi*np.fft.fftfreq(ny,d=ly/ny)
  kz = 2.*np.pi*np.fft.fftfreq(nz,d=lz/nz)
  
  kxgrid, kygrid, kzgrid, ktotsq = np.empty((nx,ny,nz)),np.empty((nx,ny,nz)),np.empty((nx,ny,nz)),np.empty((nx,ny,nz))

  for ix in range(nx):
      for iy in range(ny):
          for iz in range(nz):
              kxgrid[ix,iy,iz] = kx[ix]
              kygrid[ix,iy,iz] = ky[iy]
              kzgrid[ix,iy,iz] = kz[iz]

              ktotsq[ix,iy,iz] = kx[ix]**2 + ky[iy]**2 + kz[iz]**2
  kxgrid[int(nx/2),:,:] = 0
  kygrid[:,int(ny/2),:] = 0
  kzgrid[:,:,int(nz/2)] = 0

  ktotsq[0,0,0] = 1

  vx_k = ((H_0 * kxgrid * 1j)/ktotsq) * mass_k
  vy_k = ((H_0 * kygrid * 1j)/ktotsq) * mass_k
  vz_k = ((H_0 * kzgrid * 1j)/ktotsq) * mass_k

  vx = np.fft.ifftn(vx_k)
  vy = np.fft.ifftn(vy_k)
  vz = np.fft.ifftn(vz_k)

  return vx, vy, vz, datgrid

vx, vy, vz, smooth_r = dosmooth(3.0, Over_Den, edges[0,1]-edges[0,0], edges[1,1]-edges[1,0], edges[2,1]-edges[2,0], 1)

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

x, y, z = np.meshgrid(midx, midy, midz)
#Plot the density field and velocity field
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
#ax.scatter(x_c,y_c,z_c, alpha = 0.5)
#ax.quiver(x_c,y_c,z_c, vx_orig, vy_orig, vz_orig, alpha = 0.3, length = 0.025)
ax.quiver(x,y,z, vx, vy, vz, colors = 'k', alpha = 0.5, length = 0.01)


#X_grid = edges[0][:n_bins]
#Y_grid = edges[1][:n_bins]
#print(X_grid, Y_grid)
X,Y = np.meshgrid(edges[0][:-1],edges[1][:-1])
#X, Y = np.meshgrid(X_grid, Y_grid)

layer = np.arange(0, n_bins, 1)

for ct in layer:
    cs = ax.contourf(X,Y,smooth_r[ct], zdir = 'z', offset = edges[2][ct], cmap = plt.cm.RdYlBu_r, alpha = 0.35)
plt.colorbar(cs)
ax.set_xlabel('x-coordinate (Mpc/h)')
ax.set_ylabel('y-coordinate (Mpc/h)')
ax.set_zlabel('z-coordinate (Mpc/h)')
cs.ax.set_label('Overdensity (M_0/Mpc^3)')
ax.set_xlim(min(x_c), max(x_c))
ax.set_ylim(min(y_c), max(y_c))
ax.set_zlim(min(z_c), max(z_c))

plt.show()


