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
Simple_mass = df.Simple_Mass
Rad = df.R

#Constants
c = 3*10**5 #km/s
H_0 = 70 #km/s/Mpc
L_r0 = 0.75
L_m0 = 0.25
L_k0 = 0.044
L_l0 = 1
a = 1
f = 1

"""
data_array = np.array([[x_c], [y_c], [z_c]]).T

data_array_M = np.array([x_c,y_c,z_c, Mass_gal]).T

data_array_1 = np.array([x_c, y_c, z_c]).T
"""

#Create new mass array

x_coord = np.arange(-1000,1001, 100)
y_coord = np.arange(-1000,1001, 100)
z_coord = np.arange(-1000,1001, 100)

n_bins = len(x_coord)
print(n_bins)
#cart_coord = np.array([x_coord, y_coord, z_coord]).T
#cart_coord_mesh =np.array(np.meshgrid(x_coord, y_coord, z_coord))
#print(cart_coord_mesh[:,2,2,2])

#Mass = np.zeros((len(x_coord), len(y_coord), len(z_coord)))

#for i,x in enumerate(x_coord):
#    print(i,x)

Over_Den = np.zeros([len(x_coord),len(y_coord),len(z_coord)])
for i in range(len(x_coord)):
    for j in range(len(y_coord)):
        for k in range(len(z_coord)):
            
            if  -10<x_coord[i]<10 and -10<y_coord[j]<10 and -10<z_coord[k]<10:#np.sqrt(x_coord[i]**2 + y_coord[j]**2 + z_coord[k]**2) <= 250: : 
                Over_Den[i,j,k] = 1.0e9
            else:
                Over_Den[i,j,k] = 0#1.0e7
"""

Mass = np.zeros(len(x_coord)**3)
pos =  np.zeros((3,len(x_coord)**3))

data_array_M = np.zeros((4, len(x_coord)**3))
i = 0
for x in x_coord:
    for y in y_coord:
        for z in z_coord:
            if  np.sqrt(x**2 + y**2 + z**2) < 500: #0<x<100 and 0<y<100 and z<100:
                pos[:,i] = [x,y,z]
                Mass[i] = 10**12
                
                data_array_M[:,i] = [x,y,z, 1.0e12]
                #print(data_array_M[:,i])
            else:
                pos[:,i] = [x,y,z]
                Mass[i] = 10**9

                data_array_M[:,i] = [x,y,z,1.0e9]
                #print(data_array_M[:,i])
            i = i+1

#data_array_M = np.array([x_coord,y_coord,z_coord,Mass]).T
#print(data_array_M[:,3], data_array_M.shape)

#Create 3D histogram
#hist, edges = np.histogramdd(data_array_1, bins = (n_bins, n_bins, n_bins))
#hist, edges = np.histogramdd(pos.T, bins = (124, 124, 124))
"""

#Use these edges to define data_array_M structure

#Create density field
"""
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
    rho = den/D_Tot#(den - D_Tot) / D_Tot
    Over_Den = np.append(Over_Den, rho)


Over_Den = np.reshape(Over_Den, [(len(edges[1])-1),(len(edges[1])-1),(len(edges[1])-1)])
print(Over_Den, Over_Den.shape)
"""
"""
new = len(x_coord)
Over_Den = data_array_M[3]
Over_Den = np.reshape(Over_Den, [new, new, new])
"""


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
  

  #print('DAT',datgrid, datgrid.shape)
  
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

  vx_k = -((H_0 * kxgrid * 1j)/ktotsq) * mass_k
  vy_k = -((H_0 * kygrid * 1j)/ktotsq) * mass_k
  vz_k = -((H_0 * kzgrid * 1j)/ktotsq) * mass_k

  vx = np.fft.ifftn(vx_k)
  vy = np.fft.ifftn(vy_k)
  vz = np.fft.ifftn(vz_k)

  return vx, vy, vz, datgrid

vx, vy, vz, smooth_r = dosmooth(5.0, Over_Den, x_coord[1]-x_coord[0], y_coord[1]-y_coord[0], z_coord[1]-z_coord[0], 1)


#Mid point of bins
"""
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
"""
x, y, z = np.meshgrid(np.array([x_coord]), np.array([y_coord]), np.array([z_coord]))

#Plot the density field and velocity field
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

#ax.quiver(x,y,z, vx, vy, vx, colors = 'k', alpha = 0.5, length = 1.0e-9)

X,Y = np.meshgrid(x_coord,y_coord)
layer = np.arange(0, len(x_coord)-1, 3)
print(layer)
print(z_coord)
for ct in layer:
    print(z_coord[ct])
    cs = ax.contourf(X,Y, smooth_r[ct], zdir = 'z', offset =z_coord[ct],
                      cmap = plt.cm.RdYlBu_r, alpha = 0.3)
plt.colorbar(cs)

ax.set_xlabel('x-coordinate')
ax.set_ylabel('y-coordinate')
ax.set_zlabel('z-coordinate')
cs.ax.set_label('Overdensity')
ax.set_xlim(min(x_coord), max(x_coord))
ax.set_ylim(min(y_coord), max(y_coord))
ax.set_zlim(min(z_coord), max(z_coord))

plt.show()


