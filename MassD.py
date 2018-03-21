import numpy as np
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import math

#coordinate system assuming mass is at origin
x = ([1,2,3,4,5,6,7,8,9,10])
y = ([1,2,3,4,5,6,7,8,9,10])
print(x,y)
#biasing factor and constants
Omega = 3
f = Omega**0.55;
b = 2
rho = 0.5
i = 1
n = len(x)
#beta = f/b
beta = 1
#radial distance to points
r = np.array([])
for i in range(n):
    k = (x[i]**2 + y[i]**2)**(1/2)
    r = np.append(r,k)
    i = i+1
#print(r)

#Angle to points
theta = np.array([])
for i in range(n):
    angle = np.arcsin(r[i]/x[i]);
    theta = np.append(theta, angle);
    i = i+1
#print(theta)

#Standard radius and comparison radius
ro = 25
rs = 500

#Smoothing function with radius and peculiar velocity
W = np.array([])
for i in range(n):
    if abs(r[i]-ro)>=rs:
        Wi = 1;
    else:
        Wi = (abs(r[i]-ro)**3)/rs**3;
    W = np.append(W, Wi)
    i = i+1
print(W)

#Make all points the same weight
w = 1

inter = ([])
for i in range(n):
    intsum = sum(W[i]*w*abs(r[i]-ro)*(r[i]-ro)/(r[i]-ro)**3)+ro/3
    inter = np.append(inter, intsum)
    i = i+1
#print(inter)

u = np.array([])
for i in range(len(inter)):
    u  = beta*(((1/(4*3.14*rho))*intsum[i])+r/3);
    u = np.append(u,a);
    i = i+1
print(u)

#Components in x and y direction for the peculiar velocity
ux = np.array([])
uy = np.array([])
for i in range(n):
    ux1 = u[i]*np.cos(theta[i])
    uy1 = u[i]*np.sin(theta[i])
    ux = np.append(ux,ux1)
    uy = np.append(uy,uy1)
    i = i+1
#print(ux,uy)


#Plotting the vector field
plt.quiver(x,y,ux,uy)

plt.show()
