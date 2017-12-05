#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 15:45:30 2017


@author: robert
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def float_eq(a,b, tol=1e-12):
    return abs(a-b) < tol


def rotate(v_in, alpha, beta, gamma):
    
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    rotate_x = np.array([[1, 0, 0],
                         [0, ca, -sa],
                         [0, sa, ca]])
    
    cb = np.cos(beta)
    sb = np.sin(beta)
    rotate_y = np.array([[cb, 0, sb],
                         [0, 1, 0],
                         [-sb, 0, cb]])
    
    cg = np.cos(gamma)
    sg = np.sin(gamma)
    rotate_z = np.array([[cg, -sg, 0],
                         [sg, cg, 0],
                         [0, 0, 1]])
    
    return np.linalg.multi_dot([rotate_x, rotate_y, rotate_z, v_in])
    
def translateAndRotate(base_pts, alpha, beta, gamma, poffset):
    """
    Rotate and translate a given set of 3d coordinates
    """
    new_pts = []

    for pt in base_pts:
        new_pts.append(rotate(pt, alpha, beta, gamma) + poffset)
    
    return np.array(new_pts)   


global_cpts_list = []
current_index = 0

def appendToGlobalControlPointList(pt_array):
    
    global current_index
    global global_cpts_list
    
    connectivity = []
    new_global_pts = []
    
    # loop over all points in the user provided list coorespondong to a patgch
    for p in pt_array:
        
        pexists = False
        
        # loop over the global set of control points
        for idx, val in enumerate(global_cpts_list):
            
            # determine if the point already exists in the global array
            if float_eq(val[0],p[0]) and float_eq(val[1],p[1]) and float_eq(val[2],p[2]):
                connectivity.append(idx)
                pexists = True
                break
        
        # we get here if this is a new point not currently in the global list
        if pexists == False:
            new_global_pts.append([p[0], p[1], p[2]])
            connectivity.append(current_index)
            current_index += 1
            
    for p in new_global_pts:
        global_cpts_list.append(p)        
      
    return connectivity
    
d = 1.0             # cube dimension
t = d / 2.0        # cell thickness

print("Generating micro cell geometry with dimension %.4f  and thickness %.4f" % (d,t))

face0_pts = np.array([[0.0, 0.0, 0.0],
                    [d/2.0, 0.0, 0.0],
                    [d, 0.0, 0.0],
                    [0.0, d/2.0, 0.0],
                    [d/2.0, d/2.0, 0.0],
                    [d, d/2.0, 0.0],
                    [0.0, d, 0.0],
                    [d/2.0, d, 0.0],
                    [d, d, 0.0],
                    [t/4.0, t/4.0, t/4.0],
                    [d/2.0, t/4.0, t/4.0],
                    [d-t/4.0, t/4.0, t/4.0],
                    [t/4.0, d/2.0, t/4.0],
                    [d/2.0, d/2.0, t/4.0],
                    [d-t/4.0, d/2.0, t/4.0],
                    [t/4.0, d-t/4.0, t/4.0],
                    [d/2.0, d-t/4.0, t/4.0],
                    [d-t/4.0, d-t/4.0, t/4.0],
                    [t/2.0, t/2.0, t/2.0],
                    [d/2.0, t/2.0, t/2.0],
                    [d-t/2.0, t/2.0, t/2.0],
                    [t/2.0, d/2.0, t/2.0],
                    [d/2.0, d/2.0, t/2.0],
                    [d-t/2.0, d/2.0, t/2.0],
                    [t/2.0, d-t/2.0, t/2.0],
                    [d/2.0, d-t/2.0, t/2.0],
                    [d-t/2.0, d-t/2.0, t/2.0]])

#print("Face0 points....\n")
#for pt in face0_pts:
#    print("%.5f %.5f %.5f 1" % (pt[0],pt[1],pt[2]))
    
conn0 = appendToGlobalControlPointList(face0_pts)
print(conn0)
    
face1_pts = translateAndRotate(face0_pts, -np.pi/2.0,0.0,  0.0, [0,0,1.0])

#print("Face1 points....\n")
#for pt in face1_pts:
#    print("%.5f %.5f %.5f 1" % (pt[0],pt[1],pt[2]))

conn1 = appendToGlobalControlPointList(face1_pts)
print(conn1)
    
face2_pts = translateAndRotate(face0_pts, np.pi/2.0,0.0,  0.0, [0,1,0.0])

#print("Face2 points....\n")
#for pt in face2_pts:
#    print("%.5f %.5f %.5f 1" % (pt[0],pt[1],pt[2]))

conn2 = appendToGlobalControlPointList(face2_pts)
print(conn2)
    
face3_pts = translateAndRotate(face0_pts, 0.0, np.pi/2.0, 0.0, [0,0,1.0])

#print("Face3 points....\n")
#for pt in face3_pts:
#    print("%.5f %.5f %.5f 1" % (pt[0],pt[1],pt[2]))

conn3 = appendToGlobalControlPointList(face3_pts)
print(conn3)
    
face4_pts = translateAndRotate(face0_pts, 0.0, -np.pi/2.0, 0.0, [1,0,0.0])

#print("Face4 points....\n")
#for pt in face4_pts:
#    print("%.5f %.5f %.5f 1" % (pt[0],pt[1],pt[2]))
    
conn4 = appendToGlobalControlPointList(face4_pts)
print(conn4)

face5_pts = translateAndRotate(face0_pts, np.pi, 0.0, 0.0, [0,1,1.0])

#print("Face5 points....\n")
#for pt in face5_pts:
#    print("%.5f %.5f %.5f 1" % (pt[0],pt[1],pt[2]))

conn5 = appendToGlobalControlPointList(face5_pts)
print(conn5)


# Produce scatter plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(face0_pts[:,0],face0_pts[:,1],face0_pts[:,2])
ax.scatter(face1_pts[:,0], face1_pts[:,1], face1_pts[:,2])
ax.scatter(face2_pts[:,0], face2_pts[:,1], face2_pts[:,2])
ax.scatter(face3_pts[:,0], face3_pts[:,1], face3_pts[:,2])
ax.scatter(face4_pts[:,0], face4_pts[:,1], face4_pts[:,2])
ax.scatter(face5_pts[:,0], face5_pts[:,1], face5_pts[:,2])
plt.xlabel('x')
plt.ylabel('y')
plt.show()

conn_list = [conn0, conn1, conn2, conn3, conn4, conn5]

f = open('cube-generator-output.out','w')

f.write("3\n%d\n%d\n" % (len(global_cpts_list), len(conn_list)))

for p in global_cpts_list:
    print("%.5f %.5f %.5f 1\n" % (p[0],p[1],p[2]))
    f.write("%.5f %.5f %.5f 1\n" % (p[0],p[1],p[2]))
    
for conn in conn_list:
    f.write("2 2 2\n")
    f.write("3 3 3\n")
    f.write("0 0 0 1 1 1\n")
    f.write("0 0 0 1 1 1\n")
    f.write("0 0 0 1 1 1\n")
    for i in conn:
        f.write("%d\n" % i)     

f.close()
