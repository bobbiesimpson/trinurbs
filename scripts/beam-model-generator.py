#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:08:22 2017

@author: Robert
"""

inc = 14.0/8.0

for k in range(0,3):
    z = float(k)    
    for j in range(0,9):
        y = j * inc
        for i in range(0,3):
            x = float(i)
            print("%.10f %.10f %.10f 1" % (x,y,z))
            
knot_inc = 1.0/7.0

for i in range(0,7):
    print("%.10f " % (float(i) * knot_inc))
            
