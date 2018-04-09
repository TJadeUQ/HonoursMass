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

#Simplify to single mass and a single vector
#The mass of the objects w_i and their location r_i
w_i = ([10000000,100000000000000])
r_i = np.matrix([[100,100],[350,500]])
#Random vectors to test the field
r = np.matrix([[150,150],[150,400],[175,175],[150,100],[200,100],[50,100],[100,300],[300,300],[350,350],[200,300],[400,400],[450,400],[300,100],[250,250],[50,50]])

print(r)
print(len(r_i))
print(len(r))

#Smoothing function (Equation 4)
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
    
Smooth = SmthFunc(r_i, r, 350)
Smooth = np.matrix(Smooth)

beta = 0.5
rho = 0.008

#Equation 3- just the internal sum component using the result of the smoothing function
output = np.array([])
#Reshape the functions into arrays
Smooth = np.squeeze(np.asarray(Smooth))
r_i = np.squeeze(np.asarray(r_i))
r = np.squeeze(np.asarray(r))
print('r',r,np.shape(r))
print('r_i',r_i,np.shape(r_i))
print('Smooth',Smooth, np.shape(Smooth))

for i in range(len(r_i)):
    for j in range(len(r)):
        R = Smooth[j]*w_i[i]*((r_i[i]-r[j])/(la.norm(r_i[i]-r[j])**3))
        R = sum(R)
        output = np.append(output, R)
rshape = np.shape(r)
r_ishape = np.shape(r_i)
output = output.reshape((rshape[0],rshape[1]))

print('Ouput',output)
print(len(output))

#Equation 3 completed
Sum = np.array([])
for i in range(len(output)):
        Inter = beta*((1/(4*math.pi*rho)*(output[i]))+r[i]/3)
        print(Inter)
        Sum = np.append(Sum,Inter)

print('Sum',Sum)
#Reshape the arrays into useable forms
Sum = Sum.reshape(len(output),2)
print(np.shape(Sum))
print(Sum)
X = r[:,0]
Y = r[:,1]
print('X',X)
print(Sum[:,1])

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

