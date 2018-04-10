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
w_i = np.array([100,100,100])
r_i = np.matrix([[0,0],[500,500],[200,500]])
#Random vectors to test the field
r = np.matrix([[50,50],[100,50],[400,400],[350,200],[300,300],[200,500],[200,200],[150,100],[-500,-500]])

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

    
Smooth = SmthFunc(r_i, r, 100)
Smooth = np.matrix(Smooth)

beta = 0.5
rho = 0.008
print(np.zeros((2,len(r))))
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
        for wi,ri in zip(W,Ri):
            print(ri,wi,R[i]) #How to tell python that one is a vector and the mass is an integer
            Int = Smth[i]*wi*(ri-R[i])/(la.norm(ri-R[i])**3)          
            print(Int)
            out = np.append(out,Int)
            print(out)
    return(out)

output = Internal(r_i,r,Smooth,w_i)


print('NEW',output, output.shape)

rshape = np.shape(r)
r_ishape = np.shape(r_i)

output = np.reshape(output, (len(r),len(w_i),2))
print('Ouput',output)
print('Final',output[1,0])
print(len(output))
Final_Mix = ([])
for i,j,k in output:
    print('List',i,j,k)
    fin = i + j +k
    print('fin',fin)
    Final_Mix = np.append(Final_Mix,fin)

print('FINAL_MIX', Final_Mix, Final_Mix.shape)
Final_Mix = np.reshape(Final_Mix, (len(r),2,-1))
print(Final_Mix[0,:])
#Equation 3 completed
Sum = np.array([])
for i in range(len(Final_Mix)):
    #Inter = beta*((1/(4*math.pi*rho)*(Final_Mix[i,:])+r[i]/3))
    Inter = beta*((1/(4*math.pi*rho)*(Final_Mix[i,:])))
    print('Inter',Inter,r[i]/3,Final_Mix[i,:])
    Sum = np.append(Sum,Inter)
#The masses are forming the columns, rather than the x and y coordinates
print('Sum',Sum)
#Reshape the arrays into useable forms
Sum = Sum.reshape(len(Final_Mix),2)
print(len(Sum))

X = r[:,0]
Y = r[:,1]
print('X',X)
print('SUM',Sum)


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

