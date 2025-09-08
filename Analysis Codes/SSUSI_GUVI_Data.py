"""
Edition Date: 2025-September-01
@author: Nawapat Kaweeyanun
"""

"""
Objective: Extract data from SSUSI a netCDF file.

Prerequisite:
    SSUSI files must be saved in FolderName/dataN/Spacecraft/apl/edr-aur/year/day-of-year structure.
    GUVI files must be saved in FolderName/edr-aur/year/day-of-year structure.
"""

import netCDF4 as nc
import numpy as np
import os
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import datetime as dt

"""
Define SSUSI/GUVI extraction class.

Note: SSUSI/GUVI netCDF4 files are formatted exactly the same.
"""

class SGExtract(object):
    
    def __init__(self,filepath): #One file input at a time
        
        #Initiate netCDF4 file instance from filepath
        self.filepath = filepath
        self.nData = nc.Dataset(self.filepath)
        
        #Obtain global information from the .nc file  
        self.Instr = self.nData.MISSION #Spacecraft/Instrument e.g., 'F16'
        self.ProductType = self.nData.DATA_PRODUCT_TYPE #e.g., SDR or EDR
        self.ScanType = self.nData.SCAN_TYPE #Type of scan e.g., DISK or LIMB
        self.StartTime = self.nData.STARTING_TIME #Start time or orbit
        self.StopTime = self.nData.STOPPING_TIME #End time of orbit
        self.StartOrbNo = self.nData.STARTING_ORBIT_NUMBER #Starting orbit number
        self.StopOrbNo = self.nData.STOPPING_ORBIT_NUMBER #Final orbit number
        
        #Convert start-end times of orbit into appropriate datetimes
        self.start_dt = dt.datetime.strptime(self.StartTime,'%Y%j%H%M%S')
        self.stop_dt = dt.datetime.strptime(self.StopTime,'%Y%j%H%M%S')
        
        #Obtain variable dictionary    
        vardict = self.nData.variables
        
        #Activate data extractor methods for different product types
        if 'SDR' in self.ProductType:
            if self.ScanType == 'DISK':
                self.SDR_DISK_Extract(vardict)
            elif self.ScanType == 'LIMB':
                self.SDR_LIMB_Extract(vardict)
            else:
                print('Invalid SDR/SDR2 Mode (not Disk or Limb)')
                return
        elif 'EDR' in self.ProductType:
            if 'AURORA' in self.ProductType:
                self.EDR_Aurora_Extract(vardict)
            elif 'AURORA' not in self.ProductType and 'DISK' in self.ScanType:
                self.EDR_DISK_Extract(vardict,self.filepath)
            elif 'AURORA' not in self.ProductType and 'LIMB' in self.ScanType:
                self.EDR_LIMB_Extract(vardict,self.filepath)
            else:
                print('Data format not supported for this file.')
                return
        else:
            raise 'File not contain SDR or EDR data. No extraction.'
        
        
    def SDR_DISK_Extract(self,vardict):
        #Method for extracting SDR (Sensor Data Record) DISK data.
        #A leftover from old code. SDR data are not used for this project.
        
        #Set new dictionaries for day/night/day-auroral data
        self.SDR_Disk_day_data = {} 
        self.SDR_Disk_day_GAIM_data = {}
        self.SDR_Disk_dayauroral_data = {}
        self.SDR_Disk_dayauroral_GAIM_data = {}
        self.SDR_Disk_night_data = {} 
        self.SDR_Disk_night_GAIM_data = {}
        self.SDR_Disk_other_data = {} #Calibration/Photometer data
        
        #Looping through the variable dictionary, update correct dictionaries based on variable name filters.
        for key,value in vardict.items():
            if 'DAY' in key and 'AURORAL' not in key and 'GAIM' not in key:
                self.SDR_Disk_day_data.update({key:value})
            elif 'Day' in key and 'AURORAL' not in key and 'GAIM' in key:
                self.SDR_Disk_day_GAIM_data.update({key:value})
            elif 'DAY' in key and 'AURORAL' in key and 'GAIM' not in key:
                self.SDR_Disk_dayauroral_data.update({key:value})
            elif 'Day' in key and 'AURORAL' in key and 'GAIM' in key:
                self.SDR_Disk_dayauroral_GAIM_data.update({key:value})
            elif 'NIGHT' in key and 'GAIM' not in key:
                self.SDR_Disk_night_data.update({key:value})
            elif 'NIGHT' in key and 'GAIM' in key:
                self.SDR_Disk_night_GAIM_data.update({key:value})
            else:
                self.SDR_Disk_other_data.update({key:value})
            
        
    def SDR_LIMB_Extract(self,vardict):
        #Method for extracting SDR (Sensor Data Record) LIMB  data.
        #Same principles as the SDR DISK method above.
        #A leftover from old code. SDR data are not used for this project.
        
        self.SDR_Limb_data = {}
        self.SDR_Limb_GAIM_data = {}
        
        for key,value in vardict.items():
            if 'GAIM' not in key:
                self.SDR_Limb_data.update({key,value})
            elif 'GAIM' in key:
                self.SDR_Limb_GAIM_data.update({key,value})
                

    def EDR_Aurora_Extract(self,vardict):
        #Method for extraction Environment Data Record (EDR) Aurora data.
        #This is essentially a dummy method to rename the variable dictionary.
        #Keep variable dictionary the same as the EDR-Aurora plotting method (below) is operable based on the original data structure.
        
        self.EDR_Aurora_Data = vardict
        
    def EDR_DISK_Extract(self,vardict,filepath):
        #Method for extraction Environment Data Record (EDR) DISK data.
        #Keep variable dictionary the same, but determine if the file belongs to Day or Night modes.
        
        if 'DAY' in filepath:
            self.EDR_DayNight_String = 'DAY'
            self.EDR_Day_Disk_Data = vardict
            self.EDR_Night_Disk_Data = None
        elif 'NIGHT' in filepath:
            self.EDR_DayNight_String = 'NIGHT'
            self.EDR_Day_Disk_Data = None
            self.EDR_Night_Disk_Data = vardict

         
    def EDR_LIMB_Extract(self,vardict,filepath):
        #Method for extraction Environment Data Record (EDR) LIMB data.
        #Keep variable dictionary the same, but determine if the file belongs to Day or Night modes.
        
        if 'DAY' in filepath:
            self.EDR_DayNight_String = 'DAY'
            self.EDR_Day_Disk_Data = vardict
            self.EDR_Night_Disk_Data = None
        elif 'NIGHT' in filepath:
            self.EDR_DayNight_String = 'NIGHT'
            self.EDR_Day_Disk_Data = None
            self.EDR_Night_Disk_Data = vardict

    def EDR_Aurora_Plot(self,datadict,Chan,mode='light',save=False,savefolder=None,data_record=False):
        #Method for plotting EDR-Aurora data in AACGM projection
        #This method is not called as part of initialisation, therefore must be called separately.
        #Assume the magnetic latitude, longitude, and magnetic local times (MLT) are already in AACGM.
        
        #Obtain latitude, longitude, and disk intensity for northern and southern auroral scans.
        MagLat = datadict['LATITUDE_GEOMAGNETIC_GRID_MAP'][:,:] #Latitude
        MLT = datadict['MLT_GRID_MAP'][:,:] #Magnetic local time (same for northern and southern)
        UTN = datadict['UT_N'][:,:] #Time of scan points in UT (northern)
        UTS = datadict['UT_S'][:,:] #Time of scan points in UT (sorthern)
        
        #Generate combined intensity across selected channel(s) for northern and southern scans.
        #Set up final combined intensity array.
        IntN = np.zeros_like(MagLat)
        IntS = np.zeros_like(MagLat)
        
        #Loop over each selected channel.
        #Channels are labelled 1-5 from shortest to longest wavelengths.
        for ch in Chan:
            #Obtain the channel's radiance data for northern/southern scans.
            IntN_ch = datadict['DISK_RADIANCEDATA_INTENSITY_NORTH'][ch-1,:,:] #Convert to Python indexing
            IntS_ch = datadict['DISK_RADIANCEDATA_INTENSITY_SOUTH'][ch-1,:,:]
            
            #Set minimum threshold to filter out background.
            #This is for plotting purpose when multiple channels are combined.
            thres_min  = 100 #final minimum threshold for plotting, specifically case of combined channels
            
            #Determine specific threshold for each channels. 
            #Not every channel is calibrated the same so thresholds are different
            #If threshold is less than the minimum, update the minimum.
            if ch == 5: 
                thres = 100
                if thres < thres_min:
                    thres_min = thres 
            elif ch == 4:
                thres = 100
                if thres < thres_min:
                    thres_min = thres
            elif ch == 3:
                thres = 100
                if thres < thres_min:
                    thres_min = thres 
            elif ch == 2:
                thres = 100
                if thres < thres_min:
                    thres_min = thres 
            elif ch == 1:
                thres = 1500
                if thres < thres_min:
                    thres_min = thres 
            
            #Set all points below the threshold to zero.
            IntN_ch[IntN_ch<=thres] = 0
            IntS_ch[IntS_ch<=thres] = 0
            
            #Add the channel's intensity arrays to the total intensity arrays.
            IntN += IntN_ch
            IntS += IntS_ch
        
        #Create channel name string for file save and figure plotting.
        chan_str = 'Channel {}'.format('+'.join([str(ch) for ch in Chan]))
        
        #Record data to CSV file if prompted.
        if data_record == True:
            #Define list of parameters to be recorded (latitude, MLT, intensity-north, intensity-south).
            reclist = ['MAGLat','MLT','IntN','IntS']
            
            #Record data by parameter, row-by-row.
            for para in reclist:
                csv_filename = 'SSUSI-{}_Orb{}_{}_{}_{}_{}.csv'.format(self.Instr,self.StartOrbNo,chan_str.replace(' ',''),
                                                                        self.StartTime,self.StopTime,para)
                with open(csv_filename,'w', newline = '') as csvfile:
                    csvwriter = csv.writer(csvfile)
                    if 'lat' in para.lower():
                        for i in np.arange(len(MagLat[:,0])):
                            csvwriter.writerow(MagLat[i,:])
                    elif 'mlt' in para.lower():
                        for i in np.arange(len(MLT[:,0])):
                            csvwriter.writerow(MLT[i,:])
                    elif 'intn' in para.lower():
                        for i in np.arange(len(IntN[:,0])):
                            csvwriter.writerow(IntN[i,:])
                    elif 'ints' in para.lower():
                        for i in np.arange(len(IntS[:,0])):
                            csvwriter.writerow(IntS[i,:])
                    csvfile.close()
        
        
        #For plotting, it's easier to work with flat 1D arrays 
        MagLatFlat = MagLat.ravel()
        MLTFlat = MLT.ravel()
        IntNFlat = IntN.ravel()
        IntSFlat = IntS.ravel()
        UTNFlat = np.array(UTN.ravel()) 
        UTSFlat = np.array(UTS.ravel())
        
        #Function to convert UT decimal to hour, minute, and second. 
        def UT_Convert(UTFlat): #flat list or individual
            UT_hr = np.floor(UTFlat) #each element is float. (can be off for near midnight times?)
            UT_min = np.floor((UTFlat - UT_hr)*60)
            UT_s = np.round(((UTFlat-UT_hr)*60-UT_min)*60,0)
            
            return UT_hr,UT_min,UT_s
        
        #Convert UT for northern and southern scans.
        #All six arrays should have the same length.
        UTN_hr,UTN_min,UTN_s = UT_Convert(UTNFlat)
        UTS_hr,UTS_min,UTS_s = UT_Convert(UTSFlat)
        
        #Convert datetime from UTN/UTS to dt.datetime objects (with string arrays for plotting).
        UTN_dt = []
        UTS_dt = []
        UTN_dtstr = []
        UTS_dtstr = []
        
        #Loop over UTN_seconds flat array.
        for i in np.arange(len(UTN_s)):
            
            #Sometimes the converstion function proudce 60 for minutes and seconds.
            #Have to 'add to the next digit' and replace 60 by 0.
            if UTN_s[i] == 60:
                UTN_s[i] = 0
                UTN_min[i] = UTN_min[i]+1
            if UTS_s[i] == 60:
                UTS_s[i] = 0
                UTS_min[i] = UTS_min[i]+1                
            if UTN_min[i] == 60:
                UTN_min[i] = 0
                UTN_hr[i] = UTN_hr[i]+1            
            if UTS_min[i] == 60:
                UTS_min[i] = 0
                UTS_hr[i] = UTS_hr[i]+1
            
            #Generate datetime object for UTN and UTS inputs.
            #A provision is added in case the SSUSI/GUVI orbit straddle over midnight.
            if self.start_dt.day == self.stop_dt.day:
                Ndt = dt.datetime(self.stop_dt.year,self.stop_dt.month,self.stop_dt.day,
                                  int(UTN_hr[i]),int(UTN_min[i]),int(UTN_s[i]))
                Sdt = dt.datetime(self.stop_dt.year,self.stop_dt.month,self.stop_dt.day,
                                  int(UTS_hr[i]),int(UTS_min[i]),int(UTS_s[i]))
            elif self.start_dt.day != self.stop_dt.day:
                if UTN_hr[i] > 21: #within last 2 hours of yesterday
                    Ndt = dt.datetime(self.start_dt.year,self.start_dt.month,self.start_dt.day,
                                      int(UTN_hr[i]),int(UTN_min[i]),int(UTN_s[i]))
                elif UTN_hr[i] < 3: #within first 2 hours of today
                    Ndt = dt.datetime(self.stop_dt.year,self.stop_dt.month,self.stop_dt.day,
                                      int(UTN_hr[i]),int(UTN_min[i]),int(UTN_s[i]))
                
                if UTS_hr[i] > 21:
                    Sdt = dt.datetime(self.start_dt.year,self.start_dt.month,self.start_dt.day,
                                      int(UTS_hr[i]),int(UTS_min[i]),int(UTS_s[i]))
                elif UTS_hr[i] < 3:
                    Sdt = dt.datetime(self.stop_dt.year,self.stop_dt.month,self.stop_dt.day,
                                      int(UTS_hr[i]),int(UTS_min[i]),int(UTS_s[i]))    
            
            #Append constructed dt.datetime objects and their strings onto combined arrays.
            UTN_dt.append(Ndt)
            UTS_dt.append(Sdt)
            UTN_dtstr.append(Ndt.strftime('%Y-%m-%dT%H:%M:%S')) #convert to formatted string
            UTS_dtstr.append(Sdt.strftime('%Y-%m-%dT%H:%M:%S'))
        
        #Compute northern and southern median datetimes for plotting.
        #Repeat process from above but only for median times.
        UTNmed_hr,UTNmed_min,UTNmed_s = UT_Convert(np.nanmedian(UTNFlat[UTNFlat!=0]))
        UTSmed_hr,UTSmed_min,UTSmed_s = UT_Convert(np.nanmedian(UTSFlat[UTSFlat!=0]))
        UTNmed_hr = int(UTNmed_hr)
        UTSmed_hr = int(UTSmed_hr)
        UTNmed_min = int(UTNmed_min)
        UTSmed_min = int(UTSmed_min)
        UTNmed_s = int(UTNmed_s)
        UTSmed_s = int(UTSmed_s)
        
        if UTNmed_s == 60: #in case rounded up seconds got to 60
            UTNmed_min = UTNmed_min + 1
            UTNmed_s = 0
        if UTSmed_s == 60:
            UTSmed_min = UTSmed_min + 1
            UTSmed_s = 0

        #Save median info for later (specifically footprint overlay)
        if self.start_dt.date() == self.stop_dt.date(): #If orbit does not straddle midnight.
            self.UTNmed = self.start_dt.replace(hour = UTNmed_hr, minute=UTNmed_min, second=UTNmed_s)
            self.UTSmed = self.start_dt.replace(hour = UTSmed_hr, minute=UTSmed_min, second=UTSmed_s)
        else: #Assume northward scan always before southward scan in same orbit.
            self.UTNmed = self.start_dt.replace(hour = UTNmed_hr, minute=UTNmed_min, second=UTNmed_s)
            self.UTSmed = self.stop_dt.replace(hour = UTSmed_hr, minute=UTSmed_min, second=UTSmed_s)
        
        #Create median time strings
        UTNmed_str = '{:02d}:{:02d}:{:02d}'.format(UTNmed_hr,UTNmed_min,UTNmed_s)
        UTSmed_str = '{:02d}:{:02d}:{:02d}'.format(UTSmed_hr,UTSmed_min,UTSmed_s)
        
        #Set up plot figure, each containing northern and southern scans of the orbit in polar projections.
        fig, ax = plt.subplots(1,2, subplot_kw={'projection': 'polar'},figsize=[12,6])
        plt.subplots_adjust(wspace=0.4)
        
        #Set light and dark colour mode properties.
        if mode.lower() == 'light':
            ctick = 'k'
            tchoice = False
        elif mode.lower() == 'dark':
            ctick = 'w'
            tchoice = True
        else:
            raise Exception ('Incorrect Colour Mode')
            
        #Northern scan plot sequence. In polar projection, x = angular and y = radial axes.
        Nplot = ax[0].scatter(np.pi+(MLTFlat*np.pi/12),MagLatFlat,c=IntNFlat,
                              s=1,norm=colors.LogNorm(vmin=thres_min,vmax=10000),
                              cmap='turbo',zorder=0) #Use scatter to plot data points.
        ax[0].set_theta_zero_location('N') #Set origin to be 90 deg, not 0.
        ax[0].set_rlim(90,40,5) #Define radial latitude limit.
        xticks_loc = ax[0].get_xticks().tolist() #Obtain angular axis ticks
        ax[0].xaxis.set_major_locator(mticker.FixedLocator(xticks_loc)) #Set major axis using the obtained angular axis
        ax[0].set_xticklabels(['12','15','18','21','0','3','6','9'],fontsize=14) #Label major axis values at specific MLTs. 
        ax[0].spines['polar'].set_color(ctick) #Set angular axes edge color.
        yticks_loc = ax[0].get_yticks().tolist() #Obtain radial axis ticks.
        ax[0].yaxis.set_major_locator(mticker.FixedLocator(yticks_loc)) #Set major axis using the obtained radial axis
        ax[0].set_yticklabels(['40','','60','','80',''],fontsize=14) #Label major axis values at specific latitudes. 
        ax[0].yaxis.label.set_color(ctick) #Set radial axes ticklabel color
        ax[0].tick_params(axis='both',which='both',colors=ctick) #Set tick color for radial axis
        ax[0].grid(axis='both',which='both',color=ctick) #Set grid color
        cbar0 = fig.colorbar(Nplot,ax=ax[0],label='Disk Radiance (R)',pad=0.1) #Set colorbar.
        cbar0.outline.set_color(ctick) #Set colorbar outline color.
        cbar0.ax.axes.tick_params(axis='y',which='both',colors=ctick,labelsize=12) #Set colorbar tick color.
        cbar0.ax.axes.yaxis.label.set_color(ctick) #Set colorbar label color.
        cbar0.ax.axes.yaxis.label.set_fontsize(14) #Set colorbar label fontsize.
        ax[0].set_title('N \n' + 'Median Time = {}'.format(UTNmed_str),color=ctick,fontsize=18) #Set subplot title with median time.
        
        
        #Southern scan plot sequence (see step-by-step comments above).
        Splot = ax[1].scatter(np.pi+(MLTFlat*np.pi/12),MagLatFlat,c=IntSFlat,
                              s=1, norm=colors.LogNorm(vmin=thres_min,vmax=10000),
                              cmap='turbo',zorder=0)
        ax[1].set_theta_zero_location('N')
        ax[1].set_rlim(90,40,5)
        xticks_loc = ax[1].get_xticks().tolist()
        ax[1].xaxis.set_major_locator(mticker.FixedLocator(xticks_loc))
        ax[1].set_xticklabels(['12','15','18','21','0','3','6','9'],fontsize=14)
        ax[1].spines['polar'].set_color(ctick)
        yticks_loc = ax[1].get_yticks().tolist()
        ax[1].yaxis.set_major_locator(mticker.FixedLocator(yticks_loc))
        ax[1].set_yticklabels(['-40','','-60','','-80',''],fontsize=14)
        ax[1].yaxis.label.set_color(ctick)
        ax[1].tick_params(axis='both',which='both',colors=ctick)
        ax[1].grid(axis='both',which='both',color=ctick)
        cbar1 = fig.colorbar(Splot,ax=ax[1],label='Disk Radiance (R)',pad=0.1)
        cbar1.outline.set_color(ctick)
        cbar1.ax.axes.tick_params(axis='y',which='both',colors=ctick,labelsize=12)
        cbar1.ax.axes.yaxis.label.set_color(ctick)
        cbar1.ax.axes.yaxis.label.set_fontsize(14)
        ax[1].set_title('S \n' + 'Median Time = {}'.format(UTSmed_str),color=ctick,fontsize=18)
    
        #Set figure title, if desired.
        #fig.suptitle('Disk Radiance Orbit {} {}'.format(self.StopOrbNo,chan_str) + '\n ' 
        #             + str(self.start_dt) + ' -- ' + str(self.stop_dt),color=ctick)
        
        #Save figure, if prompted.
        if save == True:
            #Create file save name.
            savename = 'DiskRad_{}-{}-{}_Orb{}_{}.png'.format(self.stop_dt.year,self.stop_dt.month,
                                                              self.stop_dt.day,self.StopOrbNo,
                                                              chan_str.replace(' ','_'))
            #Create save directory name in Savefolder/instrument/edr-aur/year/doy/channel.
            #Define day by final datetime or the orbit.
            savedir = savefolder + '/{}/edr-aur/{}/{}/{}/'.format(self.Instr.lower(),self.stop_dt.year,
                                                               self.stop_dt.timetuple().tm_yday,chan_str)
                                                               #date/orb is always determined by final time
    
            #Make directory and save.
            os.makedirs(savedir,exist_ok=True)
            savepath = savedir + savename
            plt.tight_layout() #Use tight_layout to push figure to margins.
            plt.savefig(savepath,dpi=100,format='png',transparent=tchoice)
            plt.show
            plt.close(fig)
            return fig,ax
        
        elif save == False:
            plt.tight_layout()
            plt.close(fig)
            return fig,ax


"""
Example Script for generating all SSUSU/GUVI orbit for specific spacecraft/data-type/year/day-of-year/orbit-number.
"""

"""
#Set parameters

yr_list = [2002]
day_pick = ['077']
orbit_pick = [9]  # pick nth orbit of dayday
spacecraft = ['F16']
data_type = ['edr-aur']

#Set plot mode, channel, save, and data record properties.
mode ='light'
Chan = [3,4]
save = True
data_record = False

#Define if plotting GUVI or SSUSI so correct filepath can be generated.
satellite = 'GUVI'
if satellite == 'GUVI':
    path_start = 'Datafiles_GUVI/'
    savefolder = 'GUVI_Day_AACGM'  
elif satellite == 'SSUSI':
    path_start = 'Datafiles_SSUSI/dataN/'
    savefolder = 'SSUSI_Day_AACGM'

#Loop over each input spacecraft and data type.
for sc in spacecraft:
    for data in data_type:
        for yr in yr_list:
    
            #For each year/month layer, obtain file list, filter out non-matching day/day-of-year-orbit, then proceed down the chain.
            #Once a final filelist is obtained, call extract data from each file and call the EDR_Aurora_Plot method.
            #Use try-except to skip folders with no data files.
            
            if satellite == 'GUVI':            
                yearfolderpath = path_start + '{}/{}/'.format(data.lower(),yr)
            elif satellite == 'SSUSI':    
                yearfolderpath = path_start + '{}/apl/{}/{}/'.format(sc.lower(),data.lower(),yr)
            
            if len(day_pick) == 0:    
                daylist = os.listdir(yearfolderpath) #list of all day folders in the year given
                daylist.sort()
            else: #literally same as setting daylist = day_pick but change variable
                daylist = []
                for day in day_pick:
                    daylist.append(str(day))
            
            for day in daylist:
                dayfolderpath = yearfolderpath + day + '/'
                filelist = os.listdir(dayfolderpath)
                filelist.sort()
                
                if len(orbit_pick) != 0: #if an orbit is selected
                    filelist2 = []
                    
                    for orb in orbit_pick:
                        orb_ind = orb-1 #find python index
                        filelist2.append(filelist[orb_ind])
                    
                    filelist2.sort()
                    filelist = filelist2
                          
                os.makedirs(savefolder,exist_ok=True)
    
                for file in filelist:
                    
                    try:
                        filename = file
                        filepath = dayfolderpath + filename
                        DataSet = SGExtract(filepath)
                        fig,ax = DataSet.EDR_Aurora_Plot(DataSet.EDR_Aurora_Data,Chan,
                                                         mode,save,savefolder,data_record)
                        plt.show(fig)
                    except Exception as e: #if there is error with file opening or plotting, to next file
                        print(str(e)) 
                        print('Cannot open or plot file {} in day {} of year {}'.format(filename,day,yr))
                        
                    
                print('Plots saved for {}-{}'.format(yr,day))
"""