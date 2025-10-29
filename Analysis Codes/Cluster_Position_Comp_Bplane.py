# -*- coding: utf-8 -*-
"""
Edition Date: 2025-October-20
@author: Nawapat Kaweeyanun
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
N = np.array([87.7480,-12.3073,-137.759]) #C1 (14:54:02)

#assume SC1 is on plane, find projected position of SC2 on the plane.
RE = 6371.2 #km
SC1 = np.array([1.46293,0.12125,7.30401]) #C1 (14:54:02)
SC2 = np.array([1.44780,0.12228,7.28792]) #C3 (14:54:02)
rel_km = (SC2-SC1)*RE
plane_proj = plane_projection(SC1,SC2,N)
print(plane_proj)

