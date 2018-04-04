"""""""""""""""""""""""""""""""""
Tyler Philp
University of Queensland
43946561
Honours Project
"""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import math
#Staying in the positive quadrant
#Location of the masses
mr = ((100,200),(300,100))

k = ([10000, 20000, 30000, 40000])
#Distance to the masses ri
ri = np.array([])
for i in range(len(mr[0])):
    for j in range(len(mr[1])):
        rad = ((mr[0][i])**2+(mr[j][1])**2)**(1/2)
        ri = np.append(ri,rad)
print(k,ri)
print(np.shape(k), np.shape(ri))
#Gives 4 radii but that's fine for now

#Array of distances 
x,y = np.meshgrid(np.arange(0,400,10), np.arange(0,400,10))
r = ((x**2+y**2)**(1/2)).reshape(len(x)**2)

#Constants
rs = 100
rho = 8
beta = 1
#Peculiar Velocity function
def PecVel(Ri, R, w):
    PecVel = np.array([])
    for j in range(len(R)):
        Inter = np.array([])
        for i in range(len(Ri)):
            if abs(Ri[i]-R[j])<rs:
                SmFun = abs(Ri[i]-R[j])**3/rs**3
            else:
                SmFun = 1
            print(SmFun)
            #SmthFun = np.append(SmthFun,SmFun)

            #Incorporate summation
            Int = SmFun*(w[i]*(Ri[i]-R[j])/abs(Ri[i]-R[j])**3)
            Inter = np.append(Inter, Int)
        Sum = sum(Inter)


        #Finally PecVel!
        PV = beta*((1/(4*math.pi*rho))*(Sum)+R[j]/3)
        PecVel = np.append(PecVel, PV)
    return(PecVel)

#Test
Final = PecVel(ri,r,k)
print(Final)
print(np.shape(r))
print(np.shape(ri))

#Now need x and y components of the peculiar velocities
#Since we're in the positive quadrant, can just use cos and sin?
"""
plt.quiver(x,y,Final,-Final)
axes = plt.gca()
axes.set_aspect(1)
plt.show()
"""
