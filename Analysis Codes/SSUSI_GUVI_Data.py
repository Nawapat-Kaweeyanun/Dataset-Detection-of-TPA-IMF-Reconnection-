"""
Edition Date: 2025-October-21
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

    def EDR_Aurora_Plot(self,datadict,Chan,mode='light',cap='B',cbar_on=True,
                        save=False,savefolder=None,data_record=False):
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
                if UTN_hr[i] > 21: #Within last 2 hours of yesterday
                    if int(UTN_hr[i]) == 24: #Measure to account for orbit that started the day before but north-cap crossing occurs after midnight.
                        UTN_hr[i] = 0
                    Ndt = dt.datetime(self.start_dt.year,self.start_dt.month,self.stop_dt.day,
                                      int(UTN_hr[i]),int(UTN_min[i]),int(UTN_s[i]))
                    Sdt = dt.datetime(self.start_dt.year,self.start_dt.month,self.stop_dt.day,
                                      int(UTS_hr[i]),int(UTN_min[i]),int(UTN_s[i]))
                elif UTN_hr[i] < 3: #within first 2 hours of today
                    Ndt = dt.datetime(self.stop_dt.year,self.stop_dt.month,self.stop_dt.day,
                                      int(UTN_hr[i]),int(UTN_min[i]),int(UTN_s[i]))
                
                if UTS_hr[i] > 21:
                    if int(UTS_hr[i]) == 24: #If southern cap crossing passed after midnight. Different structure from UTN_hr[i] clause above due to differing start/stop days.
                        Sdt = dt.datetime(self.start_dt.year,self.start_dt.month,self.stop_dt.day,
                                          int(0.0),int(UTS_min[i]),int(UTS_s[i]))
                    else:
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
        UTNmed_str = self.UTNmed.strftime('%H:%M:%S')
        UTSmed_str = self.UTSmed.strftime('%H:%M:%S')
        
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
            
        #Set up figure for one subplot (north or south caps) or two subplots (both).
        if cap in ['N','S']:
            fig, ax = plt.subplots(1,1, subplot_kw={'projection': 'polar'},figsize=[4,4],squeeze=False)
            ax_ind = [0]
        elif cap == 'B':
            fig, ax = plt.subplots(1,2, subplot_kw={'projection': 'polar'},figsize=[8,4],squeeze=False)
            ax_ind = [0,1]
        else:
            raise Exception('Cap must be N, S, or B.')
        
        #Create instrument name string so no need to import new input
        if self.Instr == 'TIMED':
            instr_name = 'GUVI'
        else:
            instr_name = 'SSUSI'
            
        #Loop over 1 or 2 axis indices.     
        for ind in ax_ind:
            #Set angular axis properties.
            ax[0,ind].set_theta_zero_location('N')
            ax[0,ind].set_rlim(90,40,5)
            xticks_loc = ax[0,ind].get_xticks().tolist()
            ax[0,ind].xaxis.set_major_locator(mticker.FixedLocator(xticks_loc))
            ax[0,ind].set_xticklabels(['12','15','18','21','0','3','6','9'],fontsize=12)
            ax[0,ind].spines['polar'].set_color(ctick) #set angular axis color
            
            #Set radial axis properties and plot data.
            #Different latitude labels for northern and southern caps.
            #Because two-subplots figures are always ordered north then south, can save some if-else branches.
            yticks_loc = ax[0,ind].get_yticks().tolist()
            ax[0,ind].yaxis.label.set_color(ctick) #set y-axis ticklabel color
            ax[0,ind].tick_params(axis='both',which='both',colors=ctick) #set tick color both axes
            ax[0,ind].grid(axis='both',which='both',color=ctick) #set grid color
            if ind == 0:
                if cap in ['N','B']:
                    ax[0,ind].yaxis.set_major_locator(mticker.FixedLocator(yticks_loc))
                    ax[0,ind].set_yticklabels(['','50','60','70','',''],fontsize=10)
                    ax[0,ind].set_title('{} North\n{}'.format(instr_name,UTNmed_str) ,color=ctick,fontsize=12)
                    im = ax[0,ind].scatter(np.pi+(MLTFlat*np.pi/12),MagLatFlat,c=IntNFlat,
                                           s=1,norm=colors.LogNorm(vmin=thres_min,vmax=10000),
                                           cmap='turbo',zorder=0)
                elif cap == 'S':
                    ax[0,ind].yaxis.set_major_locator(mticker.FixedLocator(yticks_loc))
                    ax[0,ind].set_yticklabels(['','-50','-60','-70','',''],fontsize=10)
                    ax[0,ind].set_title('{} South\n{}'.format(instr_name,UTSmed_str),color=ctick,fontsize=12)
                    im = ax[0,ind].scatter(np.pi+(MLTFlat*np.pi/12),MagLatFlat,c=IntSFlat,
                                           s=1,norm=colors.LogNorm(vmin=thres_min,vmax=10000),
                                           cmap='turbo',zorder=0)
            elif ind == 1: #ind == 1 is always southern cap
                ax[0,ind].yaxis.set_major_locator(mticker.FixedLocator(yticks_loc))
                ax[0,ind].set_yticklabels(['','-50','-60','-70','',''],fontsize=10)
                ax[0,ind].set_title('{} South\n{}'.format(instr_name,UTSmed_str),color=ctick,fontsize=12)
                im = ax[0,ind].scatter(np.pi+(MLTFlat*np.pi/12),MagLatFlat,c=IntSFlat,
                                      s=1, norm=colors.LogNorm(vmin=thres_min,vmax=10000),
                                      cmap='turbo',zorder=0)
        #Set color bar properties and figure layouts.
        if cbar_on == True:
            if cap in ['N','S']:
                fig.subplots_adjust(left=0.10,right=0.72,bottom=0.06,top=0.88,wspace=0.2,hspace=0.5)
                cbar_ax = fig.add_axes([0.82,0.10,0.04,0.75]) #set axis location
                cb = plt.colorbar(im,cax=cbar_ax) #use southern plot, already set same colour limits for the two capscb.ax.tick_params(labelsize=24)
                cb.set_label(label='Radiance (R)',size=12,family='sans-serif',color=ctick)
                cb.ax.axes.tick_params(axis='y',which='both',colors=ctick) #set colorbar tick color
            elif cap == 'B':
                fig.subplots_adjust(left=0.05,right=0.85,bottom=0.01,top=0.90,wspace=0.2,hspace=0.5)
                cbar_ax = fig.add_axes([0.90,0.095,0.02,0.81]) #set axis location
                cb = plt.colorbar(im,cax=cbar_ax) #use southern plot, already set same colour limits for the two capscb.ax.tick_params(labelsize=24)
                cb.set_label(label='Radiance (R)',size=12,family='sans-serif',color=ctick)
                cb.ax.axes.tick_params(axis='y',which='both',colors=ctick) #set colorbar tick color
        else:
            plt.tight_layout()
        
        #Save figure, if prompted.
        if save == True:
            #Create file save name.
            savename = '{}_{}-{}-{}_Orb{}_{}_{}.png'.format(instr_name,self.stop_dt.year,self.stop_dt.month,
                                                            self.stop_dt.day,self.StopOrbNo,
                                                            chan_str.replace(' ','_'),cap)
            #Save folder date determined by orbit end datetime.
            savedir = savefolder + '/{}/edr-aur/{}/{}/{}/'.format(self.Instr.lower(),self.stop_dt.year,
                                                                  self.stop_dt.timetuple().tm_yday,chan_str)
    
            #Make directory and save.
            os.makedirs(savedir,exist_ok=True)
            savepath = savedir + savename
            plt.savefig(savepath,dpi=100,format='png',transparent=tchoice)
            plt.close(fig)
            return fig,ax
        elif save == False:
            plt.close(fig)
            return fig,ax