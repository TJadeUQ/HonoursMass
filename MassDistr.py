"""""""""""""""""""""""""""""""""
Tyler Philp
University of Queensland
43946561
Honours Project
"""""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt
import math

#Vectors
#Mass locations
k = 1

x,y = np.meshgrid(np.arange(-500,500,100), np.arange(-500,500,100))

#biasing factor and constants
Omega = 0.7
f = Omega**0.55;
b = 1
rho = 1.8*10**(-15) 
#beta = f/b
beta = 0.5
pi = math.pi

#Convert to spherical

r = (x**2+y**2)**(1/2)
print(np.shape(r))

ro = -5000
rs = 5000

def SmthFun(R):
    if abs(R-ro)<rs:
        SmFun = abs(R-ro)**3/rs**3
    else:
        SmFun = 1
   # print(SmFun)
    return[SmFun]

def Pec(SmFun, R):
    PecVel = beta*((1/(4*math.pi*rho))*(SmFun*k*((R-ro)/(abs(R-ro)**3)))+ro/3)
    #print(PecVel)
    return[PecVel]

Smooth = np.array([])
Peculiar = np.array([])

for i in range(x.shape[0]):
    for j in range(x.shape[1]):
        Sm = SmthFun(r[i,j])
        Smooth = np.append(Smooth, Sm)

        PecVel = Pec(Smooth[i], r[i,j])
        Peculiar = np.append(Peculiar, PecVel)
        #print(Smooth)
Peculiar = np.reshape(Peculiar, np.shape(r))
print(np.shape(Peculiar))
#Quiver needs to have scalars
plt.quiver(x,y,Peculiar)
axes = plt.gca()
axes.set_aspect(1)
plt.axis([-200,200,-200,200])
plt.show()
