"""""""""""""""""""""""""""""""""
This is the code I've written for a single mass at (0,0).
Unfortunately the function PecVel turns the peculiar velocity
into a vector of constants.
I suspect this is because I'm iterating over
position and the smoothing function simultaneously, but I'm not sure why
this would be a problem. I've reduced the number points
on the x and y axis, to simplify the problem a bit
- meshgrid was giving me a bit of stress, but I will be able to
modify the code once this other problem is sorted.

Tyler Philp
University of Queensland
43946561
Honours Project
"""""""""""""""""""""""""""""""""""





import numpy as np
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import math

#coordinate system assuming mass is at origin

xv = ([-7000,-5000,-3000,-2000,2000,3000,5000,7000])
yv = ([-7000,-5000,-3000,-2000,2000,3000,5000,7000])

#xv,yv = np.meshgrid(xv,yv)

#biasing factor and constants
Omega = 0.1
f = Omega**0.55;
b = 1
rho = 20
i = 0
n = len(xv)
#beta = f/b
beta = 1
#Convert to spherical
r = np.array([])
for i in range(n):
    c = (xv[i]**2 + yv[i]**2)**(1/2)
    r = np.append(r,c)
    i = i+1
#print(r)

theta = np.array([])
for i in range(n):
    angle = math.atan(yv[i]/xv[i]);
    theta = np.append(theta, angle);
    i = i+1
#print(theta)

ro = 2500
rs = 5000

k = 1
#Smoothing function with radius and peculiar velocity
#Define the functions for the peculiar velocity
def SmFun(R):
    SmFun = np.array([])
    for i in range(n):
        if abs(R[i]-ro)<rs:
            W = abs(R[i]-ro)**3/rs**3
        else:
            W = 1
        SmFun = np.append(SmFun, W)
        i = i+1
    return[SmFun]

def Pec(SmFun, R):
    PecVel = np.array([])
    for i in range(n):
        PV = beta*((1/(4*3.14*rho))*(SmFun[i]*k*((R[i]-ro)/(R[i]-ro)**3))+ro/3)
        PecVel = np.append(PecVel, PV)
        i = i+1
    return[PecVel]


print(r)
Smooth = SmFun(r)
Smooth = np.concatenate(Smooth)
print(Smooth)
u = Pec(Smooth, r)
u = np.concatenate(u)
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
print(ux,uy)

#Plotting the vector field
plt.quiver(xv,yv,ux,uy)

plt.show()
