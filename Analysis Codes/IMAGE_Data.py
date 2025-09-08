"""
Edition Date: 2024-August-09
@author: Nawapat Kaweeyanun
"""

"""
Objective: Extract IMAGE auroral image data from CDF file (downloadable from Cluster Science Archive)

Prerequisite: cdflib module
"""

import cdflib
import numpy as np




class data(object):
    
    def __init__(self, cdffilename, DataID):
        
        #open file with cdflib
        cdffile = cdflib.CDF(cdffilename)
        
        #Extract parameters
        Epochdates = cdffile.varget('time_tags__' + DataID)
        self.datetimes = np.array(cdflib.cdfepoch.to_datetime(Epochdates))
        
        self.hemis = cdffile.varget('Hemisphere__' + DataID) #0=North, 1=South
        if type(self.hemis) == np.int32: #remove corrupted data
            self.hemis = np.array([np.nan])
        
        self.qual = cdffile.varget('Quality__' + DataID) #1=good, 0=potential inaccuracy
        
        self.pixelnum = cdffile.varget('Pixels__' + DataID)
        
        self.image = np.array(cdffile.varget('Image__' + DataID)) #10x40x40 XY array
        
        self.XD1 = np.array(cdffile.varget('dimension_x__' + DataID)) #40-size vector, in km
        
        self.YD2 = np.array(cdffile.varget('dimension_y__' + DataID)) #40-size vector, in km
                        
