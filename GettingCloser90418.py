"""""""""""""""""""""""""""""""""
Tyler Philp
University of Queensland
43946561
Honours Project

Notes: The mass is not being iterated over correctly
Need to fix so that each one is individually calculated
"""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import linalg as la
import itertools as it
#Simplify to single mass and a single vector
#The mass of the objects w_i and their location r_i
w_i = np.array([100000000000000000000000000000000000000,10000000000])
r_i = np.matrix([[100,100],[200,500]])
#Random vectors to test the field
r = np.matrix([[150,150],[150,400],[300,450],[500,100],[400,100]])

print(r)
print(len(r_i))
print(len(r))

#Smoothing function (Equation 4)

def SmthFunc(Ri, R, rs):
    Smth = np.array([])
    for i in range(len(R)):
        for j in range(len(Ri)):
            if la.norm(Ri[j]-R[i])<rs:
                Sm_Func = (la.norm(Ri[j]-R[i])**3)/rs**3
                
            else:
                Sm_Func = 1
            Smth = np.append(Smth, Sm_Func)
    return(Smth)

    
Smooth = SmthFunc(r_i, r, 350)
Smooth = np.matrix(Smooth)

beta = 0.5
rho = 0.008

#Equation 3- just the internal sum component using the result of the smoothing function

#Reshape the functions into arrays
Smooth = np.squeeze(np.asarray(Smooth))
r_i = np.squeeze(np.asarray(r_i))
r = np.squeeze(np.asarray(r))
print('r',r,np.shape(r))
print('r_i',r_i,np.shape(r_i))
print('Smooth',Smooth, np.shape(Smooth))

def Internal(Ri,R,Smth,W):
    out = np.array([])
    for i in range(len(R)):
        for j,k in zip(W,Ri):
            print(R[i])
            print(W[j],Ri[k],R[i]) #How to tell python that one is a vector and the mass is an integer
            Int = Smth[j]*W[j]*(Ri[k]-R[i])/(la.norm(Ri[k]-R[i])**3)
            print(Int)
            Int = sum(Int)
            
            print(Int)
            out = np.append(out,Int)
    return(out)
#ERROR: only integers, slices (`:`), ellipsis (`...`),
#numpy.newaxis (`None`) and integer or boolean arrays are valid indices
output = Internal(r_i,r,Smooth,w_i)          
#The mass is not being alternated by the for loop, need to sum up the vectors  
#Want to concatenate the first two, second two and so on
print('NEW',output, output.shape)

rshape = np.shape(r)
r_ishape = np.shape(r_i)

output = np.reshape(output, (len(w_i)*rshape[0],2))
print('Ouput',output)
print(len(output))


#Equation 3 completed
Sum = np.array([])
for i in range(len(output)):
        Inter = beta*((1/(4*math.pi*rho)*(output[i]))+r[i]/3)
        print(Inter)
        Sum = np.append(Sum,Inter)
#The masses are forming the columns, rather than the x and y coordinates
print('Sum',Sum)
#Reshape the arrays into useable forms
Sum = Sum.reshape(len(output),2)
print(len(Sum))

X = r[:,0]
Y = r[:,1]
print('X',X)
print('SUM',Sum[:,0],Sum[:,1])


#Plot them vectors

fig = plt.figure()
plt.quiver(X,Y,Sum[:,0],Sum[:,1])
axes = plt.gca()
axes.set_aspect(1)
plt.ylim(0,500)
plt.xlim(0,500)
fig.suptitle('Peculiar Velocity Field')
plt.xlabel('x coordinate')
plt.ylabel('y coordinate')

plt.show()

