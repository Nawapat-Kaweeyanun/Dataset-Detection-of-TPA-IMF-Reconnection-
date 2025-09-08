"""
Edition Date: 2025-September-02
@author: Nawapat Kaweeyanun
"""

"""
Objective: Script to obtain IMAGE aurora images overlaid by Cluster footprints

Prerequisite: Cluster footprint files saved in FolderName/DataID/year/month structure
"""

import numpy as np
import datetime as dt
import SearchPlotter as SP
import OverlayPlotter as OP
from T96Analyze import T96Analyze #import class specifically
import SSUSI_GUVI_Data as SGD
import IMAGE_Timeseries as IMT


"""
Class for getting overlaying footprint data
"""

class FTOverlay(object):
    
    def __init__(self,dtimelist,sclist,FTFolder):
        self.dtimelist = dtimelist #List of datetimes
        self.sclist = sclist #List of spacecraft
        self.FTFolder = FTFolder #Path to footprint data
        
        #Create dictionary of data at selected datetimes and spacecraft
        self.dict_creator(self.dtimelist,self.sclist,self.FTFolder)

    def FT_FileExtract(self,dtime,sc,FTFolder):
        #Extract data for footprint at specified datetime and spacecraft
        FT_FileFolder = '{}/{}_CP_FGM_SPIN/{}/{}/'.format(FTFolder,sc.upper(),dtime.year,dtime.month)
        dtime_str = dtime.strftime('%Y-%m-%d_%H-%M-%S')
        FT_Filepath = FT_FileFolder + 'ClusterFootprintData__{}.csv'.format(dtime_str)
        
        #Call T96Analyze instance
        Filedata = T96Analyze()
        
        #Extract if available
        try:
            Filedata.csv_extract(FT_Filepath)
        except Exception as e:
            print(e)
            print('Footprint filepath {} not extractable'.format(FT_Filepath))
            data_list = None
            return data_list
        
        #Should only get here if try is successful
        Filedata.centralise() #Calculate mean of footprints
        Filedata.aacgm_convert() #Convert footprint to AACGM coordinates
        Filedata.mag_convert() #Convert footprint to geomagnetic coordinates (just in case)
        
        #Create data kust
        data_list = [dtime, Filedata.agm_latN, Filedata.agm_lonN, Filedata.agm_mltN, Filedata.rho_n,
                     Filedata.agm_latS, Filedata.agm_lonS, Filedata.agm_mltS, Filedata.rho_s,
                     Filedata.mag_latN, Filedata.mag_mltN, Filedata.mag_latS, Filedata.mag_mltS]
        
        print('Footprint data obtained for {} at {}'.format(sc,dtime))
        
        return data_list
    
    
    def dict_creator(self,dtimelist,sclist,FTFolder):
        #Create footprint data dictionary for given dtime list and sclist
        #Later can be broken down further by SSUSI/GUVI callers
        
        self.dt_dict = {}
        self.agm_latN_dict = {}
        self.agm_lonN_dict = {}
        self.agm_mltN_dict = {}
        self.rhoN_dict = {}
        self.mag_latN_dict = {}
        self.mag_mltN_dict = {}
        self.agm_latS_dict = {}
        self.agm_lonS_dict = {}
        self.agm_mltS_dict = {}
        self.rhoS_dict = {}
        self.mag_latS_dict = {}
        self.mag_mltS_dict = {}
        
        for sc in sclist:
            #Create empty data_list to start
            self.dt_dict.update({sc:[]})
            self.agm_latN_dict.update({sc:[]})
            self.agm_lonN_dict.update({sc:[]})
            self.agm_mltN_dict.update({sc:[]})
            self.rhoN_dict.update({sc:[]})
            self.agm_latS_dict.update({sc:[]})
            self.agm_lonS_dict.update({sc:[]})
            self.agm_mltS_dict.update({sc:[]})
            self.rhoS_dict.update({sc:[]})
            self.mag_latN_dict.update({sc:[]})
            self.mag_mltN_dict.update({sc:[]})
            self.mag_latS_dict.update({sc:[]})
            self.mag_mltS_dict.update({sc:[]})
            
            for dtime in dtimelist:
                
                #Extract data
                data_list = self.FT_FileExtract(dtime,sc,FTFolder)
                if data_list == None: #Skip if no data
                    continue
                
                #If there is data, append
                self.dt_dict[sc].append(data_list[0])
                self.agm_latN_dict[sc].append(data_list[1])
                self.agm_lonN_dict[sc].append(data_list[2])
                self.agm_mltN_dict[sc].append(data_list[3])
                self.rhoN_dict[sc].append(data_list[4])
                self.agm_latS_dict[sc].append(data_list[5])
                self.agm_lonS_dict[sc].append(data_list[6])
                self.agm_mltS_dict[sc].append(data_list[7])
                self.rhoS_dict[sc].append(data_list[8])
                self.mag_latN_dict[sc].append(data_list[9])
                self.mag_mltN_dict[sc].append(data_list[10])
                self.mag_latS_dict[sc].append(data_list[11])
                self.mag_mltS_dict[sc].append(data_list[12])
    
            #Then make array 
            self.dt_dict[sc] = np.array(self.dt_dict[sc])
            self.agm_latN_dict[sc] = np.array(self.agm_latN_dict[sc])
            self.agm_lonN_dict[sc] = np.array(self.agm_lonN_dict[sc])
            self.agm_mltN_dict[sc] = np.array(self.agm_mltN_dict[sc])
            self.rhoN_dict[sc] = np.array(self.rhoN_dict[sc])
            self.agm_latS_dict[sc] = np.array(self.agm_latS_dict[sc])
            self.agm_lonS_dict[sc] = np.array(self.agm_lonS_dict[sc])
            self.agm_mltS_dict[sc] = np.array(self.agm_mltS_dict[sc])
            self.rhoS_dict[sc] = np.array(self.rhoS_dict[sc])
            self.mag_latN_dict[sc] = np.array(self.mag_latN_dict[sc])
            self.mag_mltN_dict[sc] = np.array(self.mag_mltN_dict[sc])
            self.mag_latS_dict[sc] = np.array(self.mag_latS_dict[sc])
            self.mag_mltS_dict[sc] = np.array(self.mag_mltS_dict[sc])
            
    def dict_dtime_filter(self,dict_all,dtime):
        #Filter dictionary for datetime one wishes to plot
        plot_dict = {}
        for sc in self.sclist:
            plot_dict.update({sc:dict_all[sc][self.dt_dict[sc]==dtime]})
        return plot_dict

    def IMOverlayCaller(self,savefolder,IMtypelist=['WIC','S12'],final_save=True,mode='light',
                        saveformat='png',data_record=False,IMdt_manual_list = []):
        #Plot Cluster footprint overlaying IMAGE aurora figures.
        
        #Create union of datetimes for all spacecrafts (some do not have footprint data)
        dt_union = self.dt_dict[self.sclist[0]] #Start from one spacecraft
        for m in np.arange(1,len(self.sclist)):
            dt_union = np.union1d(dt_union,self.dt_dict[self.sclist[m]])
            
        #If union does not increase array length, then all four SCs have identical datetimes.
        if len(dt_union) == len(self.dt_dict[self.sclist[0]]):
            print('All selected spacecraft share identical datetimes.')
        else:
            print('Certain datetimes have partial spacecraft coverage.')
            
        #Separate image datetime list from footprint datetime union list, but use the same list if there is no manual input
        if len(IMdt_manual_list) == 0:
            IMdt_union = dt_union.copy()
        else:
            IMdt_union = IMdt_manual_list
                
        #Impose max image limit, based on letter availability
        if len(IMdt_union) > 26:
            raise Exception('More than 26 subplots requested. Stop before labelling error.')
        
        #design number of subplot columns by divisor
        if len(IMdt_union) == 1:
            ncol = 1
        elif len(IMdt_union) <= 4:
            ncol = len(dt_union)
        elif len(IMdt_union) in [5,13,17,21,25]:
            #Special cases where 3 looks best
            ncol = 3
        else:
            ncol = 4
        
        for IMtype in IMtypelist:
            #Create IMAGE plot
            IMfig,IMax,hno_arr = IMT.timeseries_plot(IMdt_union, ncol, IMtype,None,False,mode,'png',False) #Not save intially.
            #hno is number 0 = north, 1 = south, np.nan = else/corrupted
            
            #Call overlay function using footprint dictionaries as inputs
            IMfig,IMax = OP.IMTSOverlay(IMfig,IMax,IMtype,hno_arr,dt_union,self.dt_dict,self.agm_latN_dict,
                                        self.agm_mltN_dict,self.rhoN_dict,self.agm_latS_dict,self.agm_mltS_dict,
                                        self.rhoS_dict,final_save,savefolder,mode,saveformat,data_record)
            
            
    
    def SGOverlayCaller(self,satellite,datafolder,savefolder,SSinstrlist=['F16','F17','F18','F19'],
                        Chan=[5],final_save=True,mode='light',saveformat='png',
                        filepath_manual_list = []):
        
        #Plot Cluster footprint overlaying SSUSI/GUVI aurora figures.
        
        #Get all file paths within footprint datetime range (only if no manual input)
        if len(filepath_manual_list) == 0: 
            filepath_list = []
            
            #Search for closest-orbit SSUSI/GUVI data file
            for dtime in self.dtimelist:
                if satellite == 'SSUSI':
                    for SSinstr in SSinstrlist:
                        filepath = SP.SSUSI_filesearch(datafolder,'EDR-Aur',SSinstr,dtime)
                        filepath_list.append(filepath)
                elif satellite == 'GUVI':
                    filepath = SP.GUVI_filesearch(datafolder, 'EDR-Aur', dtime)
                    filepath_list.append(filepath)
        
            #Remove file duplicates
            filepath_list = list(set(filepath_list))
        #Or can input data file path
        else:
            filepath_list = filepath_manual_list
        
        
        #Extract SSUSI/GUVI data then plot, then overlay
        for filepath in filepath_list:
            
            #Extract data
            DataSet = SGD.SGExtract(filepath)
            #Plot aurora scan
            fig,ax = DataSet.EDR_Aurora_Plot(DataSet.EDR_Aurora_Data,Chan,mode,save=False)
            
            #For each footprint datetime, overlay footprint
            for ft_dt in self.dtimelist:
                ft_latN = self.dict_dtime_filter(self.agm_latN_dict,ft_dt) #change from agm to mag
                ft_mltN = self.dict_dtime_filter(self.agm_mltN_dict,ft_dt)
                ft_rhoN = self.dict_dtime_filter(self.rhoN_dict,ft_dt)
                ft_latS = self.dict_dtime_filter(self.agm_latS_dict,ft_dt)
                ft_mltS = self.dict_dtime_filter(self.agm_mltS_dict,ft_dt)
                ft_rhoS = self.dict_dtime_filter(self.rhoS_dict,ft_dt)
                    
                #Note: toggle switch == True allows footprint to be flipped (see switch function in OverlayPlotter.py)
                if satellite == 'SSUSI': 
                    for SSinstr in SSinstrlist:
                        fig,ax = OP.SGOverlay(fig,ax,ft_dt,ft_latN,ft_mltN,ft_rhoN,ft_latS,ft_mltS,
                                              ft_rhoS,satellite,SSinstr,Chan,DataSet.StopOrbNo,
                                              True,mode,final_save,savefolder,saveformat)
                elif satellite == 'GUVI':
                    fig,ax = OP.SGOverlay(fig,ax,ft_dt,ft_latN,ft_mltN,ft_rhoN,ft_latS,ft_mltS,
                                          ft_rhoS,satellite,'TIMED',Chan,DataSet.StopOrbNo,
                                          True,mode,final_save,savefolder,saveformat)

"""
Main Script

Uncomment IMAGE or SSUSI/GUVI section depending on use

Footprint data directory is the name shown in the directory this file is located.

"""
"""

sclist = ['C1']
FTFolder = '<FT Data Directory Name>'
dtimelist = [dt.datetime(2002,3,18,14,15,0),
             dt.datetime(2002,3,18,14,35,0),
             dt.datetime(2002,3,18,14,55,0),
             dt.datetime(2002,3,18,15,15,0)]
FT = FTOverlay(dtimelist,sclist,FTFolder)
"""
"""
satellite = 'IMAGE'
IMtypelist = ['WIC']
save = True
savefolder='IMAGE_Timeseries_T96'
mode = 'light'
saveformat = 'png'
data_record =False
IMdt_manual_list = []
FT.IMOverlayCaller(savefolder,IMtypelist,save,mode,saveformat,data_record,
                   IMdt_manual_list)
"""

"""
satellite = 'GUVI'
datafolder = 'Datafiles_{}'.format(satellite)
SSinstrlist = ['F16']
Chan = [3,4]
save = True
savefolder = '{}_T96'.format(satellite)
mode = 'light'
saveformat = 'png'
filepath_manual_list = ['Datafiles_GUVI/edr-aur/2002/113/']
FT.SGOverlayCaller(satellite,datafolder,savefolder,SSinstrlist,Chan,save,mode,saveformat,
                   filepath_manual_list)
"""
