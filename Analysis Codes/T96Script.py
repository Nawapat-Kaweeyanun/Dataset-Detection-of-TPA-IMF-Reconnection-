"""
Edition Date: 2025-September-01
@author: Nawapat Kaweeyanun
"""
"""
Objective: From a comma-limited file containing Cluster magnetometer (FGM) data,
operate Tsyganenko-96 (T96) trace model to obtain Cluster's footprint location

Prerequisites:
    1) CSV Cluster files (written using CSV function in CSA_Instr_Data.py) with Cluster FGM data in spin resolution.
    2) T96 model installed in Python framework (see Coxon-Larquier (2020))

Note: Call script provided for application. See details below.    
"""

import os
import numpy as np
import tsyganenko as tsy
import csv
import re
import datetime as dt
import CDAWeb_Extract as CDAEx



"""
Class for obtaining Cluster's footprint through Tsyganenko-96 model
"""


class FootprintFinder(object):
    
    def __init__(self, csvfilepath): 
        #Class called for each CSV filepath.
        self.csvfilepath = str(csvfilepath)
        
    def getdata(self):
        #Extract data from CSV file.
        self.datetime = [] #Datetime
        self.hint = [] #Half-interval
        self.Bvec = [] #Magnetic field vector
        self.Bmag = [] #Magnetic field magnitude
        self.SC_gse = [] #Spacecraft position in GSE
        self.FGMrange = [] #Fluxgate magnetometer range
        self.FGMtele = [] #Fluxgate magnetoemter telemetry
    
        #Open csv file and read data row by row.
        with open(self.csvfilepath) as csvfile:
            dataset = csv.reader(csvfile, delimiter = ',')
            datetime_init = [] #initial datetime to be split then converted
            hint_init = []
            Bvec_init = []
            Bmag_init = []
            SC_gse_init = []
            FGMrange_init = []
            FGMtele_init = []
            
            for row in dataset:
                datetime_init.append(row[0])
                hint_init.append(row[1])
                Bvec_init.append(row[2])
                Bmag_init.append(row[3])
                SC_gse_init.append(row[4])
                FGMrange_init.append(row[5])
                FGMtele_init.append(row[6])
            
            #Remove headers
            datetime_init.pop(0)
            hint_init.pop(0)
            Bvec_init.pop(0)
            Bmag_init.pop(0)
            SC_gse_init.pop(0)
            FGMrange_init.pop(0)
            FGMtele_init.pop(0)
            
            #Convert to arrays
            datetime_init = np.array(datetime_init)
            hint_init = np.array(hint_init)
            Bvec_init = np.array(Bvec_init)
            Bmag_init = np.array(Bmag_init)
            SC_gse_init = np.array(SC_gse_init)
            FGMrange_init = np.array(FGMrange_init)
            FGMtele_init = np.array(FGMtele_init)
            
            #Convert datetime from string to dt.datetime object
            for row in np.arange(len(datetime_init)):
                dt_str = datetime_init[row]
                
                #Fill in missing microseconds if needed. Microseconds is preceded by a dot that non-microseconds items do not have.
                if '.' not in dt_str:
                    ('a')
                    dt_str = dt_str + '.000000'
                
                #Obtain dt.datetime object from formatted string, then append it to datetime list.
                dtime = dt.datetime.strptime(dt_str,'%Y-%m-%d %H:%M:%S.%f')
                self.datetime.append(dtime)
            
            #Convert datetime list to array.
            self.datetime = np.array(self.datetime)
            
            #Convert other parameters from strings to float
            for row in np.arange(len(hint_init)):
                #Half interval
                self.hint.append(float(hint_init[row]))
                
                #For magnetic field vector, must strip brackets and extra spaces and split the components first.
                Bvec_edit = Bvec_init[row].replace('[','').replace(']','') #remove brackets
                Bvec_edit2 = re.sub(' +',' ',Bvec_edit).strip(' ') 
                Bvec_split = Bvec_edit2.split(' ')
                self.Bvec.append([float(Bvec_split[0]),float(Bvec_split[1]),
                                 float(Bvec_split[2])])
                
                #Magnetic field magntitude
                self.Bmag.append(float(Bmag_init[row]))
                
                #For spacecraft position, process data as magnetic field vector above.
                SC_edit = SC_gse_init[row].replace('[','').replace(']','')
                SC_edit2 = re.sub(' +',' ',SC_edit).strip(' ')
                SC_split = SC_edit2.split(' ')
                self.SC_gse.append([float(SC_split[0]),float(SC_split[1]),
                                    float(SC_split[2])])
                
                #FGM range
                self.FGMrange.append(FGMrange_init[row])
                
                #FGM telemetry
                self.FGMtele.append(FGMtele_init[row])

            #Convert lists to arrays.
            self.hint = np.array(self.hint)
            self.Bvec = np.array(self.Bvec)
            self.Bmag = np.array(self.Bmag)
            self.SC_gse = np.array(self.SC_gse)
            self.FGMrange = np.array(self.FGMrange)
            self.FGMtele = np.array(self.FGMtele)
    
    def para_opt(self, vsw_gse=[-400., 0., 0.], Pdyn=2, dst=-5, By_imf=0,
                 Bz_imf=-5, l_max=5000, rmax=60, rmin=1, dsmax=0.01, 
                 err=0.000001):
        #Set optional parameters for T96 trace class (see Coxon Python wrapper).
        self.vsw_gse = vsw_gse #Solar wind velocity in GSE
        self.Pdyn = Pdyn #Dynamic pressure
        self.dst = dst #DST index
        self.By_imf = By_imf #IMF By
        self.Bz_imf = Bz_imf #IMF Bz
        self.l_max = l_max #Maximum number of tracing steps
        self.rmax = rmax #Maximum radius for tracing
        self.rmin = rmin #Minimum radius for tracing
        self.dsmax = dsmax #Maximum step size
        self.err = err #Tracing step tolerance.
        
    def startcoord_convert(self):
        #Method to convert starting coordinate necessary for T96 model input.
        #See Geopack Recalc_08 and Coxon Python wrapper.
        
        #Preliminary check to make sure that array lengths are equal. (So can make loop).
        if len(self.datetime) != len(self.SC_gse[:,0]):
            raise "Datetime array length different from coordinate array length"
        
        self.rho_geo = np.zeros_like(self.SC_gse[:,0]) #In km
        self.lat_geo = np.zeros_like(self.SC_gse[:,0])
        self.lon_geo = np.zeros_like(self.SC_gse[:,0])
        
        #Loop over datetime.
        for ip in np.arange(len(self.datetime)):
            
            #Call the recalculation function.
            tsy.geopack.recalc_08(self.datetime[ip].year,
                                  self.datetime[ip].timetuple().tm_yday,
                                  self.datetime[ip].hour,
                                  self.datetime[ip].minute,
                                  self.datetime[ip].second, *self.vsw_gse)
            
            #Convert coordinates from GSE to GSW then to geocentric (GEO), then to GEO polar coordinates.
            xgsw,ygsw,zgsw = tsy.convert.coordinates(self.SC_gse[ip,0],
                                                     self.SC_gse[ip,1],
                                                     self.SC_gse[ip,2],
                                                     "GSE","GSW")
            xgeo,ygeo,zgeo = tsy.convert.coordinates(xgsw,ygsw,zgsw,"GSW","GEO")
            rgeo,colatgeo,longeo = tsy.convert.car_to_sph(xgeo,ygeo,zgeo)
            latgeo = 90-(colatgeo*180/np.pi)
            longeo = longeo*180/np.pi
            self.rho_geo[ip] = rgeo
            self.lat_geo[ip] = latgeo
            self.lon_geo[ip] = longeo
 
    def getfootprint(self,ft_datetime,savedir):
        #Method to generate footprint for a list of trace datetime lists. ('snapshots')
        
        #Loop over the trace datetime list.
        for row in np.arange(len(ft_datetime)):

            #Check if the trace datetime is among the datetimes in the footprint file 
            #Count a match if it is within two seconds of trace datetime (should produce only one match for FGM 4-s spin resolution) 
            delta_s_thres = 2;
            
            #Create list to append all matches
            window_datetime = []
            window_index = []
            
            #Loop over the footprint file's datetime list and append all datetimes that match the trace datetime.
            for i in np.arange(len(self.datetime)):
                delta = abs(ft_datetime[row] - self.datetime[i])
                if delta.days == 0 and delta.seconds < delta_s_thres:
                    window_datetime.append(self.datetime[i])
                    window_index.append(True)                    
                else:
                    window_index.append(False)
                    
                
            #Convert lists to arrays. 
            window_datetime = np.array(window_datetime)
            window_index = np.array(window_index)
            
            #Filter spacecraft coordinates (polar geocentric) for datetimes to be traced.
            window_lat_geo = self.lat_geo[window_index]
            window_lon_geo = self.lon_geo[window_index]
            window_rho_geo = self.rho_geo[window_index]
            
            #Skip file if no match is found (determine using latitude array)
            if len(window_lat_geo) == 0:
                print('No available data in window at ft_datetime no. ' + str(row+1))
                continue
            
            Coord = "GEO"
            
            #If there is only one datetime, remove list.
            if len(window_datetime) == 1: #special case?
                window_datetime = window_datetime[0]
            
            #Call Trace function to obtain footprint.
            WindowTrace = tsy.Trace(window_lat_geo,window_lon_geo,window_rho_geo,
                                    Coord,window_datetime,self.vsw_gse,
                                    self.Pdyn,self.dst,self.By_imf,self.Bz_imf,
                                    self.l_max,self.rmax,self.rmin,self.dsmax,
                                    self.err)
            
            print('Footprint Traced')

            #Write footprint location in geocentric polar coordinates to a csv file.
            snapname = str(ft_datetime[row]).replace(':','-').replace('.','-')
            snapname = snapname.replace(' ','_')
            csvname = 'ClusterFootprintData__' + snapname + '.csv'
            
            #Define headers for csv file
            fields = ['lat_n','lon_n ','rho_n','lat_s','lon_s','rho_s']
            
            #Make save directory and construct filepath to save.
            os.makedirs(savedir,exist_ok=True)
            csvpath = savedir + csvname
            
            #Open file to write data
            with open(csvpath, 'w', newline = '') as csvfile:
                csvwr = csv.writer(csvfile)
                
                csvwr.writerow(fields)
                
                
                for n in np.arange(len(WindowTrace.lat_n)):
                    csvwr.writerow([WindowTrace.lat_n[n],WindowTrace.lon_n[n],
                                    WindowTrace.rho_n[n],WindowTrace.lat_s[n],
                                    WindowTrace.lon_s[n],WindowTrace.rho_s[n]])
                    
            print('CSV file saved at {}'.format(savedir))
        


"""
Example Script for T96 Footprint Tracing
"""
#Set initial parameters
sc_choice = ['C1']
PathStart = 'Datafiles_ClusterFGM_CSV'
SaveTarget = 'Cluster_FT_T96_New'

#Set datetime range and interval
dti = dt.datetime(2002,3,18,14,15,0)
dtf = dt.datetime(2002,3,18,15,16,0)
interval = 20 #in minutes

#Create footprint datetime list using while loop
dt_iter = dti
ft_datetime = [dt_iter]
while dt_iter < dtf:
    dt_iter = dt_iter + dt.timedelta(seconds=60*interval)
    ft_datetime.append(dt_iter)
 

for sc in sc_choice:
    DataID = '{}_CP_FGM_SPIN'.format(sc.upper())
        
    #Set up filepaths down to month (use start datetime in case of cross-month tracing)
    FileFolderPath = '{}/{}/{}/{}/'.format(PathStart,DataID,dti.year,dti.month)
    Filelist = os.listdir(FileFolderPath)
    Filelist.sort()
    
    #Apply day filter based to obtain file matching start date
    Filelist2 = []
    for file in Filelist:
        if '{}-{:02d}-{:02d}'.format(dti.year,dti.month,dti.day) in file.split('_')[2]:
            Filelist2.append(file)
    Filelist2.sort()
    Filelist = Filelist2 #filter down the file list.
    
    #Loop over all the files obtained.
    for file_no in np.arange(len(Filelist)):
        #Obtain file path
        Filename = Filelist[file_no]
        Filepath = FileFolderPath + '/' + Filename
        
        #Initialize the Finder class (see above) for this file
        CFT = FootprintFinder(Filepath)
        CFT.getdata()
        #CFT.para_opt()
        
        #Compute average solar wind conditions for datetime range defined by the file.
        SWdti = CFT.datetime[0] #First datetime in file
        SWdti = SWdti.strftime('%Y-%m-%dT%H:%M:%SZ') #Formatted to string
        SWdtf = CFT.datetime[-1] #Last datetime in file
        SWdtf = SWdtf.strftime('%Y-%m-%dT%H:%M:%SZ') #Formatted to string    
        paralist = ['By','Bz','Vx','Vy','Vz','P'] #Average IMF Y/Z, SW velocity, and SW pressure
        OG = 'OMNI (Combined 1AU IP Data; Magnetic and Solar Indices)' #Define OMNI grids and coordinates
        IT = 'Magnetic Fields (space)'
        res = '1min'
        grid = 'HRO2'
        coord = 'GSE' #Note: the Coxon Python wrapper should operate in GSE.
        bundle = CDAEx.CDAWeb_Manager(SWdti,SWdtf,paralist,OG,IT,res,grid,coord) #Call OMNI data from CDAWeb/
        
        #Calculate mean SW parameter values (ignoring NaNs)
        para_meandict = {}
        for para in bundle.varlist: #Note: parameter names differ from paralist above.
            bundle_para = bundle.data[para]
            para_mean = np.nanmean(bundle_para)
            para_meandict.update({para:para_mean})
                
        #Set solar wind condition inputs (unspecifieds are defaults)
        CFT.para_opt(vsw_gse=[para_meandict['Vx'],para_meandict['Vy'],para_meandict['Vz']],
                     Pdyn=para_meandict['Pressure'],By_imf=para_meandict['BY_GSE'],Bz_imf=para_meandict['BZ_GSE'])
        
        #If there is spacecraft position data in this file, skip.                       
        if len(CFT.SC_gse) == 0: 
            continue
        
        #Convert coordinates before calling T96 model.        
        CFT.startcoord_convert()
        
        #Define directory for outputs to be saved.
        SaveFolderPath = '{}/{}/{}/{}/'.format(SaveTarget,DataID,dti.year,dti.month)
        
        #Trace footprint for this file
        CFT.getfootprint(ft_datetime,SaveFolderPath)
        print('file no. ' + str(file_no+1) + ' scouted')
        