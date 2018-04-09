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

w_i = ([10000000,1000000])
r_i = np.matrix([[100,100],[500,500]])

r = np.matrix([[300,300],[350,350],[200,300],[400,400],[450,400],[300,100],[250,250],[50,50]])

print(r)
print(len(r_i))
print(len(r))

diff = np.array([])
for i in range(len(r_i)):
    for j in range(len(r)):
        di = la.norm(r_i[i]-r[j])
        print(di)
        diff = np.append(diff, di)

print(diff)
diff = np.matrix(diff)
print(diff)
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

    
Smooth = SmthFunc(r_i, r, 200)
Smooth = np.matrix(Smooth)

beta = 0.5
rho = 0.008

#Now want to find the function
#Convert this to a callable function
#As well as introduce the constants
output = np.array([])

Smooth = np.squeeze(np.asarray(Smooth))
r_i = np.squeeze(np.asarray(r_i))
r = np.squeeze(np.asarray(r))
print('r',r,np.shape(r))
print('r_i',r_i,np.shape(r_i))
print('Smooth',Smooth)

for i in range(len(r_i)):
    for j in range(len(r)):
        R = Smooth[i]*w_i[i]*((r_i[i]-r[j])/(la.norm(r_i[i]-r[j])**3))
        output = np.append(output, R)
rshape = np.shape(r)
r_ishape = np.shape(r_i)
output = output.reshape((rshape[0],r_ishape[0],rshape[1]))

print('Ouput',output)
print(len(output))

#Now try to do the sum

Sum = np.array([])
for i in range(len(output)):
        Inter = beta*((1/(4*math.pi*rho)*sum(output[i]))+r[i]/3)
        print(Inter)
        Sum = np.append(Sum,Inter)

print(Sum)

Sum = Sum.reshape(2,len(output))

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

