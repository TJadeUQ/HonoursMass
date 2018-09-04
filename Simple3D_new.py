import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
from scipy import linalg as la
import array
from mpl_toolkits import mplot3d
import pandas as pd
"""

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
Simple_mass = df.Simple_Mass
Rad = df.R
"""
#Constants
c = 3*10**5 #km/s
H_0 = 70 #km/s/Mpc
L_r0 = 0.75
L_m0 = 0.25
L_k0 = 0.044
L_l0 = 1
a = 1
f = 1

#Create new mass array

x_coord = np.arange(-1000,1001, 275)
y_coord = np.arange(-1000,1001, 275)
z_coord = np.arange(-1000,1001, 275)
#512 bins
n_bins = len(x_coord)
print(n_bins)

Over_Den = np.zeros([len(x_coord),len(y_coord),len(z_coord)])
for i in range(len(x_coord)):
    for j in range(len(y_coord)):
        for k in range(len(z_coord)):
            if np.sqrt(x_coord[i]**2 + y_coord[j]**2 + z_coord[k]**2) < 500:
                Over_Den[i,j,k] = 5.0e8
            else:
                Over_Den[i,j,k] = 1.0e8

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

vx, vy, vz, smooth_r = dosmooth(1.0, Over_Den, x_coord[1]-x_coord[0], y_coord[1]-y_coord[0], z_coord[1]-z_coord[0], 1.0)
x, y, z = np.meshgrid(np.array([x_coord]), np.array([y_coord]), np.array([z_coord]))

#Plot the density field and velocity field
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

ax.quiver(x,y,z, vx, vy, vz, colors = 'k', alpha = 0.5, length = 1.0e-9)

X,Y = np.meshgrid(2*x_coord,2*y_coord)
layer = np.arange(0, len(z_coord), 1)
for ct in layer:
    cs = ax.contourf(X,Y, smooth_r[ct], zdir = 'z', offset = z_coord[ct],
                      cmap = plt.cm.RdYlBu_r, alpha = 0.5)
plt.colorbar(cs)

ax.set_xlabel('x-coordinate')
ax.set_ylabel('y-coordinate')
ax.set_zlabel('z-coordinate')
cs.ax.set_label('Overdensity')
ax.set_xlim(2*min(x_coord),2*max(x_coord))
ax.set_ylim(2*min(y_coord),2*max(y_coord))
ax.set_zlim(2*min(z_coord),2*max(z_coord))

plt.show()


