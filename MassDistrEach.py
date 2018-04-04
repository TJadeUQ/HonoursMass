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
mr = ((10,20),(30,10))
k = ([10000, 20000, 30000, 40000])
#Distance to the masses ri
ri = np.array([])
for i in range(len(mr[0])):
    for j in range(len(mr[1])):
        rad = ((mr[0][i])**2+(mr[j][1])**2)**(1/2)
        ri = np.append(ri,rad)

#Gives 4 radii but that's fine for now

#Array of distances 
x,y = np.meshgrid(np.arange(1,1000,100), np.arange(1,1000,100))
r = ((x**2+y**2)**(1/2)).reshape(len(x)**2)
x = x.reshape(100)
y = y.reshape(100)
theta = np.array([])
for i in range(len(r)):
    angle = math.acos(x[i]/r[i])
    theta = np.append(theta,angle)
    i = i+1
print(np.shape(theta))
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

ThankyaJesus = np.array([])
for i in range(len(r)):
    Sumterm = Sum(ri,r[i],k)
    Peculiar = PecVel(Sumterm, k, r[i], ri)
    ThankyaJesus = np.append(ThankyaJesus, Peculiar)
    i = i+1
print(np.shape(ThankyaJesus))
ThankyaJesus = ThankyaJesus.reshape(4,100)
print(np.shape(ThankyaJesus))
#Now need x and y components of the peculiar velocities
#Since we're in the positive quadrant, can just use cos and sin?
Pecx = np.array([])
Pecy = np.array([])
for i in range(len(ThankyaJesus)):
    pecx = ThankyaJesus[i]*np.cos(theta[i])
    pecy = ThankyaJesus[i]*np.sin(theta[i])
    Pecx = np.append(Pecx, pecx)
    Pecy = np.append(Pecy, pecy)
    i = i+1


plt.quiver(x,y,Pecx,Pecy)
axes = plt.gca()
#axes.set_aspect(1)
plt.show()
