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
w_i = (100,200,300,400)
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
output = output.reshape((20,2))
print(output)

#measure from the origin
origin = [0],[0]
#Plot them vectors
plt.quiver(*origin, output[:,0], output[:,1])
plt.show()
