import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
import array
from mpl_toolkits import mplot3d
import pandas as pd


df = pd.read_csv('3377.csv')
z_red = df.z
x_c = df.x_c
y_c = df.y_c
z_c = df.z_c
vx_orig = df.vx
vy_orig = df.vy
vz_orig = df.vz
vr = df.vr
log_m = df.log_m
dec = df.dec
ra = df.ra
Mass_gal = df.Mass
print(max(x_c),max(y_c),max(z_c))
#Constants
c = 3*10**5 #km/s
H_0 = 70 #km/s/Mpc
L_r0 = 0.75
L_m0 = 0.25
L_k0 = 0.044
L_l0 = 1
a = 1
f = 1
n_bins = 16

#Want to subsample data
Mass_gal_sub = Mass_gal[::2]
x_c_sub = x_c[::2]
y_c_sub = y_c[::2]
z_c_sub = z_c[::2]
print(len(Mass_gal_sub))



#Create data array
def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', 0.0e0)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value
    return vector

data_array = np.array([x_c_sub, y_c_sub, z_c_sub]).T

#Create 3D histogram, zero pad data
hist, edges = np.histogramdd(data_array, bins = (n_bins,n_bins,n_bins))
edges = np.array(edges)

delta_g = (hist/np.mean(hist))-1

edges_x = np.pad(edges[0], (int(n_bins/2),int(n_bins/2)), 'reflect', reflect_type = 'odd')
edges_y = np.pad(edges[1], (int(n_bins/2),int(n_bins/2)), 'reflect', reflect_type = 'odd')
edges_z = np.pad(edges[2], (int(n_bins/2),int(n_bins/2)), 'reflect', reflect_type = 'odd')
edges = np.array([edges_x,edges_y,edges_z])

mass = np.pad(delta_g, int(n_bins/2), pad_with)

scale_x = edges[0,1]- edges[0,0]
scale_y = edges[1,1] - edges[1,0]
scale_z = edges[2,1] - edges[2,0]
scale_factor = np.average((scale_x,scale_y, scale_z))





# Smooth gridded distribution - Chris Blake code
def dosmooth(scale,datgrid,lx,ly,lz,b):

  nx,ny,nz = datgrid.shape[0],datgrid.shape[1],datgrid.shape[2]
  norm = np.sum(datgrid)
  
  x = lx*np.fft.fftfreq(nx, d = lx/nx)
  y = ly*np.fft.fftfreq(ny, d = ly/ny)
  z = lz*np.fft.fftfreq(nz, d = lz/nz)

  rsqgrid = x[:,np.newaxis,np.newaxis]**2 + y[np.newaxis,:,np.newaxis]**2 + z[np.newaxis,np.newaxis,:]**2 #Transform into radial coordinates
  kerngrid = np.sqrt(2.0 * np.pi * scale**2) * np.exp(-2.0 * np.pi**2 * rsqgrid**2 * scale**2) #Gaussian smoothing kernel
  #kerngrid = np.fft.fft(kerngrid)
  datspec = np.fft.fftn(datgrid)
  #kernspec = np.fft.fftn(kerngrid)#Fourier transform of smoothed kernel
  #mass_k = datspec*kernspec
  mass_k = datspec * kerngrid
  
  datgrid = b * np.fft.ifftn(mass_k)
  datgrid = datgrid * norm/np.sum(datgrid)

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
              
  kxgrid[int(nx/2),:,:] = 0.0
  kygrid[:,int(ny/2),:] = 0.0
  kzgrid[:,:,int(nz/2)] = 0.0

  ktotsq[0,0,0] = 1.0
    
  vx_k = -((H_0 * kxgrid * 1j)/ktotsq) * mass_k
  vy_k = -((H_0 * kygrid * 1j)/ktotsq) * mass_k
  vz_k = -((H_0 * kzgrid * 1j)/ktotsq) * mass_k

  vx = np.fft.ifftn(vx_k)
  vy = np.fft.ifftn(vy_k)
  vz = np.fft.ifftn(vz_k)

  return vx, vy, vz, datgrid

vx, vy, vz, smooth_r = dosmooth(scale_factor/150.0, mass, edges[0,2*n_bins]-edges[0,0], edges[1,2*n_bins]-edges[1,0], edges[2,2*n_bins]-edges[2,0], 1.0)


#Mid point of bins
midx = np.zeros(2*n_bins)
midy = np.zeros(2*n_bins)
midz = np.zeros(2*n_bins)
for I in range(2*n_bins):       
    midx[I] = (edges[0,I+1]+edges[0,I])/2
for J in range(2*n_bins):
    midy[J] = (edges[1,J+1]+edges[1,J])/2
for K in range(2*n_bins):
    midz[K] = (edges[2,K+1]+edges[2,K])/2
    
#exit()
x, y, z = np.meshgrid(midx, midy, midz)
#Plot the density field and velocity field
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

ax.scatter(x_c_sub,y_c_sub,z_c_sub, alpha = 0.5)
#ax.quiver(x_c,y_c,z_c, vx_orig, vy_orig, vz_orig, alpha = 0.3, length = 0.005)
#ax.quiver(x,y,z, -vx, -vy, -vz, colors = 'k', alpha = 0.5, length = 0.00005)

X,Y = np.meshgrid(edges[0][:-1],edges[1][:-1])
layer = np.arange(0, 2*n_bins, 5)

for ct in layer:
    cs = ax.contourf(X,Y,smooth_r[ct], zdir = 'z', offset = edges[2][ct], cmap = plt.cm.RdYlBu_r, alpha = 0.35)

plt.colorbar(cs)
ax.set_xlabel('x-coordinate (Mpc/h)')
ax.set_ylabel('y-coordinate (Mpc/h)')
ax.set_zlabel('z-coordinate (Mpc/h)')
cs.ax.set_label('Overdensity (M_0/Mpc^3)')
ax.set_xlim(min(edges[0]), max(edges[0]))
ax.set_ylim(min(edges[1]), max(edges[1]))
ax.set_zlim(min(edges[2]), max(edges[2]))

plt.show()
