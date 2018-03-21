import numpy as np
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import math

#coordinate system assuming mass is at origin
x = np.linspace(-50.0,50.0)
y = np.linspace(-50.0,50.0)

xv,yv = np.meshgrid(x,y)

#biasing factor and constants
Omega = 0.1
f = Omega**0.55;
b = 1
rho = .5
i = 0
n = len(x)


r = (xv**2 + yv**2)**(1/2)
print(r)

#theta = np.array([])
#for i in range(n):
 #   angle = math.atan(yv[i]/xv[i]);
  #  theta = np.append(theta, angle);
  #i = i+1
theta = np.arcsin(yv/r)
print(theta)

ro = 25
rs = 25

#Smoothing function with radius and peculiar velocity
#W = np.array([])
#for i in range(n):
#if abs(r-ro)>=rs:
 #   W = 1;
#else:
W = (abs(r-ro)**3)/rs**3;

    #W = np.append(W, Wi)
    #i = i+1
print(W)

w = 1

#inter = ([])
#for i in range(n):
intsum = sum(W*w*abs(r-ro)*(r-ro)/(r-ro)**3)+ro/3
    #inter = np.append(inter, intsum)
    #i = i+1
print(intsum)

#u = np.array([])
#for i in range(len(inter)):
u  = (f/b)*(((1/(4*3.14*rho))*intsum)+r/3);
#    u = np.append(u,a);
#    i = i+1
#print(u)

#Components in x and y direction for the peculiar velocity
#ux = np.array([])
#uy = np.array([])
#for i in range(n):
ux = u*np.cos(theta)
uy = u*np.sin(theta)
    #ux = np.append(ux,ux1)
    #uy = np.append(uy,uy1)
    #i = i+1
print(ux,uy)

#Smoothing function with redshift
#if abs(z-zo)>=zs:
 #   W = 1
#else:
 #   W = (abs(z-zo).^3)./zs^3


#Plotting the vector field
#plt.axis('equal')

plt.quiver(xv,yv,ux,uy)

plt.show()
