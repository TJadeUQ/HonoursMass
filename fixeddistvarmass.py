"""""""""""""""""""""""""""""""""
Tyler Philp
University of Queensland
43946561
Honours Project
Mass vs velocity for single distance
Notes:
"""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import linalg as la
import itertools as it
#Simplify to single mass and a single vector
#The mass of the objects w_i and their location r_i
w_i = np.array(np.linspace(0,1000,num = 10))
r_i = np.array([[250,250]])
#Random vectors to test the field
r = np.array([100,100])

def SmthFunc(Ri, R, rs, W):
    Smth = np.array([])
    for i in range(len(W)):
        for j in range(len(Ri)):
            if la.norm(Ri[j]-R)<rs:
                Sm_Func = (la.norm(Ri[j]-R)**3)/rs**3
                
            else:
                Sm_Func = 1
            Smth = np.append(Smth, Sm_Func)
    return(Smth)
Smooth = np.array(SmthFunc(r_i,r,500,w_i))
print(Smooth)
beta = 0.5
rho = 0.008

Smooth = np.squeeze(np.asarray(Smooth))
r_i = np.squeeze(np.asarray(r_i))
r = np.squeeze(np.asarray(r))
print('r',r,np.shape(r))
print('r_i',r_i,np.shape(r_i))
print('Smooth',Smooth, np.shape(Smooth))


def Internal(Ri,R,Smth,W):
    out = np.array([])
    for wi,Sm in zip(W,Smth):
        print(Ri,wi,Sm) 
        Int = Sm*wi*(Ri-R)/(la.norm(Ri-R)**3)          
        print('HERE',Int)
        out = np.append(out,Int)
        print(out)
    return(out)

output = Internal(r_i,r,Smooth,w_i)


print('NEW',output, output.shape)

rshape = np.shape(r)
r_ishape = np.shape(r_i)

output = np.reshape(output, (len(w_i),2))
print('Ouput',output)
print('Final',output[1,0])

Sum = np.array([])
for i in range(len(output)):
    #Inter = beta*((1/(4*math.pi*rho)*(Final_Mix[i,:])+r[i]/3))
    Inter = beta*((1/(4*math.pi*rho)*(output[i,:])))
    print('Inter',Inter,r/3,output[i,:])
    Sum = np.append(Sum,Inter)

print('Sum',Sum.shape)
Sum = np.reshape(Sum,(len(w_i),2))
print('Sum',Sum)
mag_pec = np.array([])
for i in Sum:
    magSum = la.norm(i)
    print(magSum)
    mag_pec = np.append(mag_pec,magSum)
print(w_i.shape,mag_pec.shape)

#Center the graph at zero by using the differenc between the mass and the vector
plt.scatter(w_i,mag_pec)
plt.title('Peculiar velocity at (100,100) for a variety of masses')
plt.xlabel('Mass')
plt.ylabel('Magnitude of peculiar velocity')

plt.show()
