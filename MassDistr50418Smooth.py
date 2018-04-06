"""""""""""""""""""""""""""""""""
Tyler Philp
University of Queensland
43946561
Honours Project

Notes: Need to include constants
Some of the vectors are looking a bit shonky
Will need to see more vectors
Want to avoid meshgrid at all costs
"""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import linalg as la
#mass and r of the ith galaxy
w_i = (1000,100)
r_i = np.array([[10.0,20.0],[20.0,10.0]])
r_ishape = np.shape(r_i)
print(r_ishape)

#The r vectors - eventually want a meshgrid? - or just arbitary vectors?
#Leaving them as vectors makes life a lot easier
x = np.arange(100)
y = np.arange(100)
r = np.array([])
for i in range(len(x)):
    coords = tuple([x[i-10],y[i]])
    coords1 = tuple([x[i-25],y[i-10]])
    coords2 = tuple([x[i-30],y[i-15]])
    coords3 = tuple([x[i-45],y[i-20]])
    r = np.append(r, coords)
    r = np.append(r, coords1)
    r = np.append(r, coords2)
    r = np.append(r, coords3)

r = r.reshape(4*len(x),2)

rshape = np.shape(r)

#Smoothing function
def SmthFunc(Ri, R, rs):
    Smth = np.array([])
    for i in range(len(Ri)):
        for j in range(len(R)):
            if la.norm(Ri[i]-R[j])<rs:
                Sm_Func = (la.norm(Ri[i]-R[j])**3)/rs**3
            else:
                Sm_Func = 1
            print(Sm_Func)
            Smth = np.append(Smth, Sm_Func)
    return(Smth)
    
Smooth = SmthFunc(r_i, r, 15)
print('Smooth', Smooth)
beta = 1
rho = 0.008
#Now want to find the function
#Convert this to a callable function
#As well as introduce the constants
output = np.array([])
for i in range(len(r_i)):
    for j in range(len(r)):
        R = beta*(1/(4*math.pi*rho)*(Smooth[i]*w_i[i]*((r_i[i]-r[j])/(la.norm(r_i[i]-r[j])**3)))+r[j]/3)
        output = np.append(output, R)
output = output.reshape((rshape[0],r_ishape[0],rshape[1]))

print(output)


#Now try to do the sum
Sum_X = np.array([])
Sum_Y = np.array([])
for i in range(len(output)):
        Inter_X = sum(output[i,0])
        Inter_Y = sum(output[i,1])
        Sum_X = np.append(Sum_X, Inter_X)
        Sum_Y = np.append(Sum_Y, Inter_Y)

print(Sum_X.shape, Sum_Y.shape)

X = r[:,0]
Y = r[:,1]
print(X.shape,Y.shape)


#Plot them vectors
fig = plt.figure()
plt.quiver(X,Y,Sum_X,Sum_Y)
axes = plt.gca()
axes.set_aspect(1)
fig.suptitle('Peculiar Velocity Field')
plt.xlabel('x coordinate')
plt.ylabel('y coordinate')

plt.show()
