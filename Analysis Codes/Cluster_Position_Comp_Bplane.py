# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 18:29:21 2024

@author: nk1c22
"""

"""
Objective: Compare two spacecraft positions in plane normal to magnetic field (reference frame) measured by the first spacecraft

"""

import numpy as np

def plane_projection(RSC1,RSC2,N):
    #projection of SC2 position on plane perpendicular to SC1
    #N must be measured at position RSC1
    RE = 6371.2
    rel_km = (RSC2-RSC1)*RE #this is SC2 position in Cartesian originated at SC1

    Nsq = np.sum(N**2) #square each element then sum up array

    vec_proj = ((np.dot(rel_km,N))/Nsq)*N #in km
    plane_proj = rel_km-vec_proj
    return plane_proj

"""

Main Script

"""


#define plane normal vectors
N = np.array([93.8400,-8.1686,-139.0670]) #C1 (14:55:02)

#assume SC1 is on plane, find projected position of SC2 on the plane.
RE = 6371.2 #km
SC1 = np.array([1.491728,0.108878,7.316323]) #C1 (14:55:02)
SC2 = np.array([1.476566,0.113125,7.300205]) #C3 (14:55:01)
rel_km = (SC2-SC1)*RE
plane_proj = plane_projection(SC1,SC2,N)
print(plane_proj)

