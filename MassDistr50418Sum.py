"""""""""""""""""""""""""""""""""
Tyler Philp
University of Queensland
43946561
Honours Project

Notes: New problems are coming down to how to define r
How to get python to take in vectors to produce a vector
The definition of the coordinate system is the problem at this point
"""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt
import math

#mass and r of the ith galaxy
w_i = (1000,200,300,400)
r_i = np.array([[1,1],[2,2],[6,6],[10,10]])
print(np.shape(r_i))

#The r vectors - eventually want a meshgrid
r = np.array([[0,0],[1,2],[10,8],[8,8],[3,4]])
print(np.shape(r))

#Now want to find ri-r for each case
output = np.array([])
for i in range(len(r_i)):
    for j in range(len(r)):
        R = w_i[i]*(r_i[i]-r[j])
        output = np.append(output, R)
output = output.reshape((5,4,2))
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
