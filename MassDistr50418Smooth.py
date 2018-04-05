"""""""""""""""""""""""""""""""""
Tyler Philp
University of Queensland
43946561
Honours Project
"""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import linalg as la
#mass and r of the ith galaxy
w_i = (1000,200,300,400)
r_i = np.array([[1,1],[2,2],[6,6],[10,10]])
r_ishape = np.shape(r_i)
print(r_ishape)

#The r vectors - eventually want a meshgrid? - or just arbitary vectors?
#Leaving them as vectors makes life a lot easier
r = np.array([[0,0],[1,5],[1,2],[10,8],[8,8],[3,4],[5,5]])
rshape = np.shape(r)
print(rshape[0])


#Smoothing function
def SmthFunc(Ri, R, rs):
    Smth = np.array([])
    for i in range(len(Ri)):
        for j in range(len(R)):
            if la.norm(Ri[i]-R[j])<rs:
                Sm_Func = (la.norm(Ri[i]-R[j])**3)/rs**3
            else:
                Sm_Func = 1
        Smth = np.append(Smth, Sm_Func)
    return(Smth)
    
Smooth = SmthFunc(r_i, r, 3)
print(Smooth)

#Now want to find ri-r for each case
output = np.array([])
for i in range(len(r_i)):
    for j in range(len(r)):
        R = Smooth[i]*w_i[i]*((r_i[i]-r[j])/(la.norm(r_i[i]-r[j])**3))+r[j]/3
        output = np.append(output, R)
output = output.reshape((rshape[0],r_ishape[0],rshape[1]))
print(output)
print(len(output))
print(output[:,0],output[:,1])
#Now try to do the sum
Sum_X = np.array([])
Sum_Y = np.array([])
for i in range(len(output)):
        Inter_X = sum(output[i,0])
        Inter_Y = sum(output[i,1])
        Sum_X = np.append(Sum_X, Inter_X)
        Sum_Y = np.append(Sum_Y, Inter_Y)
print(Sum_X,Sum_Y)
print(Sum_X.shape, Sum_Y.shape)

x = r[:,0]
y = r[:,1]
print(x.shape,y.shape)
print(x,y)
#Plot them vectors
plt.quiver(x,y,Sum_X,Sum_Y)
plt.show()
