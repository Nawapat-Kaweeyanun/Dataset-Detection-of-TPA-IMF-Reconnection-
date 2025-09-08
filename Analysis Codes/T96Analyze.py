"""
Edition Date: 2024-August-09
@author: Nawapat Kaweeyanun
"""

"""
Objective: Store information from T96-traced Cluster footprint data file.

Prerequisite: aacgmv2 module
"""

import csv
import numpy as np
from spacepy import coordinates as coords
from spacepy.time import Ticktock
from spacepy.toolbox import rad2mlt
import datetime as dt
import aacgmv2

"""
Define class for manipulating T96 data (single file)
"""

class T96Analyze(object):
    
    def __init__(self):
        self.RE = 6371.2
    
    def csv_extract(self,filepath):
        #Extract csv file from given path

        #Set up lists to fill with data
        lat_n_str = []
        lon_n_str = []
        rho_n_str = []
        lat_s_str = []
        lon_s_str = []
        rho_s_str = []
        
        #Open CSV file and append variables, one row at a time.
        with open(filepath) as csvfile:
            CData = csv.reader(csvfile,delimiter = ',')
            for row in CData:
                lat_n_str.append(row[0])        
                lon_n_str.append(row[1])
                rho_n_str.append(row[2])
                lat_s_str.append(row[3])
                lon_s_str.append(row[4])
                rho_s_str.append(row[5])

        #Remove headers
        lat_n_str.pop(0)
        lon_n_str.pop(0)
        rho_n_str.pop(0)
        lat_s_str.pop(0)
        lon_s_str.pop(0)
        rho_s_str.pop(0)

        #Set up NumPy arrays then convert strings to floats
        self.lat_n = []
        self.lon_n = []
        self.rho_n = []
        self.lat_s = []
        self.lon_s = []
        self.rho_s = []
        self.hemis = []

        for i in np.arange(len(lat_n_str)):
            self.lat_n.append(float(lat_n_str[i]))
            self.lon_n.append(float(lon_n_str[i]))
            self.rho_n.append(float(rho_n_str[i]))
            self.lat_s.append(float(lat_s_str[i]))
            self.lon_s.append(float(lon_s_str[i]))
            self.rho_s.append(float(rho_s_str[i]))
    
        self.lat_n = np.array(self.lat_n)
        self.lon_n = np.array(self.lon_n)
        self.rho_n = np.array(self.rho_n)
        self.lat_s = np.array(self.lat_s)
        self.lon_s = np.array(self.lon_s)
        self.rho_s = np.array(self.rho_s)
   
        dt_str = filepath.split('/')[-1].split('__')[1].split('.')[0]
        self.dtime = dt.datetime.strptime(dt_str,'%Y-%m-%d_%H-%M-%S')
        self.dtime_str = self.dtime.strftime('%Y-%m-%dT%H:%M:%S')
        
    def centralise(self):
        #Call to find means if csv values are effectively the same
        #e.g., all within one second interval
        self.lat_n = np.array([np.mean(self.lat_n)])
        self.lon_n = np.array([np.mean(self.lon_n)])
        self.rho_n = np.array([np.mean(self.rho_n)])
        self.lat_s = np.array([np.mean(self.lat_s)])
        self.lon_s = np.array([np.mean(self.lon_s)])
        self.rho_s = np.array([np.mean(self.rho_s)])
        
    def aacgm_convert(self):
        #Convert coordinates from geocentric to AACGM
        
        #Calculate altitude above Earth's surface for northern footprint
        altN = self.rho_n-self.RE
        for row in np.arange(len(altN)):
            if altN[row] < 0 and altN[row] > -10**-2: #remove all numerical-effect negatives
                altN[row] = 0
        
        #Convert northern footprint to AACGM coordinates
        ft_agm_N = aacgmv2.get_aacgm_coord_arr(self.lat_n,self.lon_n,altN,self.dtime,method='GEOCENTRIC')
        self.agm_latN = ft_agm_N[0]
        self.agm_lonN = ft_agm_N[1]
        self.agm_mltN = ft_agm_N[2]
        
        
        #Calculate altitude above Earth's surface for southern footprint
        altS = self.rho_s-self.RE
        for row in np.arange(len(altS)):
            if altS[row] < 0 and altS[row] > -10**-2: #remove all numerical-effect negatives
                altS[row] = 0

        #Convert southern footprint to AACGM coordinates
        ft_agm_S = aacgmv2.get_aacgm_coord_arr(self.lat_s,self.lon_s,altS,self.dtime,method='GEOCENTRIC')
        self.agm_latS = ft_agm_S[0]
        self.agm_lonS = ft_agm_S[1]
        self.agm_mltS = ft_agm_S[2]
        
        self.dtime_arr = np.array([self.dtime]*len(self.agm_latN))
        
        
    def mag_convert(self):
        #Convert coordinates from geocentric to geomagnetic at a particular datetime
        
        #Set up datetime array of the same size as the footprint arrays
        dtime_arr = [self.dtime_str]*len(self.rho_n)
        
        #Set up an instance of geocentric coordinate object in SpacePy, fixedd datetime to each footprint, then convert to geomagnetic.
        #(for northern footprint)
        ft_geoN = coords.Coords(np.transpose([self.rho_n,self.lat_n,self.lon_n]),'GEO','sph',units=['km','deg','deg'])
        ft_geoN.ticks = Ticktock(dtime_arr,'UTC')
        ft_magN = ft_geoN.convert('MAG','sph') #['km', 'deg', 'deg']
        
        #Obtain geomagnetic coordinates data including MLT
        self.mag_rhoN = ft_magN.data[:,0]
        self.mag_latN = ft_magN.data[:,1]
        self.mag_lonN = ft_magN.data[:,2]
        self.mag_mltN = rad2mlt(self.mag_lonN*np.pi/180)

        #Similar process for southern footprint as northern footprint
        ft_geoS = coords.Coords(np.transpose([self.rho_s,self.lat_s,self.lon_s]),'GEO','sph',units=['km','deg','deg'])
        ft_geoS.ticks = Ticktock(dtime_arr,'UTC')
        ft_magS = ft_geoS.convert('MAG','sph')
        
        self.mag_rhoS = ft_magS.data[:,0]
        self.mag_latS = ft_magS.data[:,1]
        self.mag_lonS = ft_magS.data[:,2]
        self.mag_mltS = rad2mlt(self.mag_lonS*np.pi/180)
        