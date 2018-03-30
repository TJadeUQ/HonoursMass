"""""""""""""""""""""""""""""""""
Tyler Philp
University of Queensland
43946561
Honours Project
"""""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import math

#Vectors

x = [-2000,2000,-2000,2000,-5000,5000,-5000,5000,-1000,-1000,-3000,-3000,-2500,0,3000,1000,-6000,6000,-3000,-5500]
y = [-2000,2000,2000,-2000,-5000,5000,5000,-5000,-1000,1000,-3000,3000,-2500,1000,2000,5000,-6000,-3000,-5000,0]

#xv,yv = np.meshgrid(xv,yv)

#biasing factor and constants
Omega = 0.7
f = Omega**0.55;
b = 1
rho = 1.8*10**(-13) 
i = 0
n = len(x)
#beta = f/b
beta = 1
pi = math.pi
#Convert to spherical
r = np.array([])
for i in range(len(x)):
    radial = -(x[i]**2+y[i]**2)**(1/2)
    r = np.append(r, radial)

    i = i+1
print(r)

theta = np.array([])
for i in range(len(x)):
    if x[i]>0 and y[i]>=0:
        angle = math.atan(y[i]/x[i])
        
    elif x[i]>=0 and y[i]<0:
        angle = 3*pi/2-math.atan(x[i]/y[i])
        
    elif x[i]<=0 and y[i]>0:
        angle = pi/2 - math.atan(x[i]/y[i])

    else:
        angle = pi + math.atan(y[i]/x[i])
        
    theta = np.append(theta, angle)
    i = i+1

    
#Constants - distance to the object ro, density smoothing radius - rs, weight - k
ro = 2500
rs = 5000

k = np.ones(n)
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
        PV = beta*((1/(4*math.pi*rho))*(SmFun[i]*k[i]*((R[i]-ro)/(abs(R[i]-ro)**3)))+ro/3)
        print(PV)
        PecVel = np.append(PecVel, PV)
        i = i+1
    return[PecVel]



Smooth = SmFun(r)
Smooth = np.concatenate(Smooth)

u = Pec(Smooth,r)
u = np.concatenate(u)


#Components in x and y direction for the peculiar velocity
ux = np.array([])
uy = np.array([])
for i in range(n):
    ux1 = u[i]*np.cos(theta[i])
    uy1 = u[i]*np.sin(theta[i])
    ux = np.append(ux,ux1)
    uy = np.append(uy,uy1)
    i = i+1


#Plotting the vector field
plt.quiver(x,y,ux,uy)
axes = plt.gca()
axes.set_aspect(1)
plt.show()
