import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import cProfile
import statistics as stat
from anytree import Node, RenderTree, NodeMixin

#Import the data
df = pd.read_csv('Order.csv')

z = df.z
x_c = df.x_c
y_c = df.y_c
z_c = df.z_c
vx = df.vx
vy = df.vy
vz = df.vz
vr = df.vr
log_m = df.log_m
Mass_gal = df.Mass
"""
Create the arrays
Want to find the medians and then
be able to attribute those medians to the
nodes - the leaves should be each each object
"""
#Create data set of redshift and mass
z = np.array(z)
Mass_gal = np.array(Mass_gal)
data = np.stack((z,Mass_gal), axis = 1)
#Create dictionary to associate values with each other
data_lib = dict(zip(z,Mass_gal))

def find_nearest(array,value):
    idx = (np.abs(array[:,0]-value)).argmin()
    return array[idx]

#Now find the median of the left and right child
#Should really makes these callable functions

#Now want to take child_left and child_right and run them as new roots
#Will have the problem that the median is no longer a value in the dictionary

#Create class of Node
class MyBaseClass(object):
    foo = 4
class MyClass(MyBaseClass, NodeMixin):
    def __init__(self, name, z_coord, Mass_coord, parent = None):
        super(MyClass, self).__init__()
        self.name = name
        self.z_coord = z_coord
        self.Mass_coord = Mass_coord
        self.parent = parent

#Create the root
z_mid = stat.median(data[:,0])
Mass_gal_mid = data_lib[z_mid]

arr = np.array([z_mid, Mass_gal_mid])

med_data = stat.median(data[:,0])
mid_data = find_nearest(data, med_data)

Root = MyClass('Root',z_mid, Mass_gal_mid)

def left_right_array(arr,data):
    left = ([])
    right = ([])
    for i in data:
        if i[0] < arr[0]:
            left  = np.append(left, [i])
        else:
            right = np.append(right, [i])
    left.shape = (-1,2)
    right.shape = (-1,2)
    return(left,right)


def children(L,R,parent):
    mid_left = find_nearest(L, stat.median(L[:,0]))
    mid_right = find_nearest(R, stat.median(R[:,0]))

    Left_Child = MyClass('Left_Child', mid_left[0], mid_left[1], parent = parent)
    Right_Child = MyClass('Right_Child', mid_right[0], mid_right[1], parent = parent)

    for pre, _, node in RenderTree(parent):
        treestr = u"%s%s" % (pre, node.name)
        print(treestr.ljust(8), node.z_coord, node.Mass_coord)

    return(Left_Child, Right_Child, parent,mid_left,mid_right)  

left_array,right_array = left_right_array(arr,data)
left_child,right_child,adult,mid_left,mid_right = children(left_array,right_array,Root)

#Want to iterate until len(data) = 1, want to iterate so I pass the next iteration the left and right arrays as the new data
#will need to iteratire over 2^n nodes with each pass

#Try the next pass (Next node on the left)
left_2, right_2 = left_right_array(mid_left, left_array)


left_child_2, right_child_2,adult_2,mid_left_2,mid_right_2 = children(left_2, right_2, left_child)


#Want to try and run this for each of the nodes in the layer 2^(depth)
"""
for i in len(2^(Left_Child.depth):
    run the code
    do the thing
"""
