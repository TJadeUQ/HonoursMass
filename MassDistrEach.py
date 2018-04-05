"""""""""""""""""""""""""""""""""
Tyler Philp
University of Queensland
43946561
Honours Project

Notes: So instead of creating four masses, this code is generating 4 vectors per point
Need to seperate the idea of ri and r
The summation term is also complicating things
Want to get to the point where I can have a field of r then input some ri and wi
to give a uniqie field - might be worthwhile to think about this in reverse
Start with the r field of uniform vectors then add the masses
"""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import math
#Staying in the positive quadrant
#Location of the masses
mr = ([100,250],[300,10])
k = ([10000, 20000])

#Distance to the masses ri
ri = np.array([])
for i in range(len(mr[0])):
    print(mr[0][i])
    rad = ((mr[0][i])**2+(mr[1][i])**2)**(1/2)
    ri = np.append(ri,rad)
print(ri)
#Gives 4 radii but that's fine for now

#Array of distances and convert to polar 
x,y = np.meshgrid(np.arange(1,1000,100), np.arange(1,1000,100))
r = ((x**2+y**2)**(1/2)).reshape(len(x)**2)
print(r)
x = x.reshape(100)
y = y.reshape(100)
theta = np.array([])
for i in range(len(r)):
    angle = math.acos(x[i]/r[i])
    theta = np.append(theta,angle)
    i = i+1
print((theta))
#x = 20
#y = 2000
#r = ([(x**2+y**2)**(1/2)])
#Constants
rs = 450
rho = 0.8
beta = 1

#Smoothing function needs to be individual for each r as well as summing over the masses - need to
#do the sum first

def Sum(Ri,R,w):
    Sum = np.array([])
    for i in range(len(Ri)):
        if abs(Ri[i]-R)<rs:
            SmFun = abs(Ri[i]-R)**3/rs**3
        else:
            SmFun = 1
        S = (SmFun*w[i]*(Ri[i] - R)/(abs(Ri[i]-R))**3)
        Sum = np.append(Sum,S)
        i = i+1
    return(Sum)

def PecVel(Sum, w, R, Ri):
    PV = np.array([])
    for i in range(len(Ri)):
        u = beta*(((1/(4*math.pi*rho))*sum(Sum))+R/3)
        PV = np.append(PV, u)
        i = i+1
    return(PV)

Final = np.array([])
for i in range(len(r)):
    Sumterm = Sum(ri,r[i],k)
    Peculiar = PecVel(Sumterm, k, r[i], ri)
    Final = np.append(Final, Peculiar)
    i = i+1
print((Final))

#Now need x and y components of the peculiar velocities
#Since we're in the positive quadrant, can just use cos and sin?
print(np.shape(Final))
print(np.shape(theta))

Pecx = np.array([])
Pecy = np.array([])
for i in range(len(theta)):
    pecx = Final[i]*np.cos(theta[i])
    pecy = Final[i]*np.sin(theta[i])
    Pecx = np.append(Pecx, pecx)
    Pecy = np.append(Pecy, pecy)
    i = i+1


plt.quiver(x,y,Pecx,Pecy)
axes = plt.gca()
#axes.set_aspect(1)
plt.show()
