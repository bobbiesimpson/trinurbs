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


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

    
def translateAndRotate(base_pts, axis, theta, poffset):
    """
    Rotate and translate a given set of 3d coordinates
    """
    new_pts = []

    for pt in base_pts:
        new_pts.append(np.dot(rotation_matrix([1,0,0], -np.pi/2.0), pt) + poffset)
    
    return np.array(new_pts)   

    
d = 1.0             # cube dimension
t = d / 10.0        # cell thickness

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

for pt in face0_pts:
    print("%.5f %.5f %.5f 1" % (pt[0],pt[1],pt[2]))
    
face1_pts = translateAndRotate(face0_pts, [1,0,0], -np.pi/2.0, [0,0,1.0])

print("Face1 points....\n")
for pt in face1_pts:
    print("%.5f %.5f %.5f 1" % (pt[0],pt[1],pt[2]))
    
face2_pts = translateAndRotate(face1_pts, [1,0,0], 2 * np.pi, [0,1,0.0])

print("Face2 points....\n")
for pt in face2_pts:
    print("%.5f %.5f %.5f 1" % (pt[0],pt[1],pt[2]))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(face0_pts[:,0],face0_pts[:,1],face0_pts[:,2])
ax.scatter(face1_pts[:,0], face1_pts[:,1], face1_pts[:,2])
ax.scatter(face2_pts[:,0], face2_pts[:,1], face2_pts[:,2])
plt.show()





