"""
Edition Date: 2025-August-29
@author: Nawapat Kaweeyanun
"""


"""
Objective: Overlay SSUSI/GUVI on top of IMAGE (greyscaled)

Rerequisite: Run from same directory as IMAGE and SSUSI/GUVI data directories

"""

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator
import SearchPlotter as SP
import IMAGE_Data as IMD
import SSUSI_GUVI_Data as SGD
import OverlayPlotter as OP
from Footprint_Aurora_Overlay import FTOverlay


def SGonIM(dtime,IMtype,SGind=0,SGChan=[5],mfl=[],mode='light',save=False,saveformat='png'): #mfl = manual file input list
    #Function to overlay SSUSI/GUVI scan on grayscaled IMAGE auroral picture.
    #ACCEPT ONLY ONE DATETIME!

    #Define data folder names.
    IMFolder = 'Datafiles_IMAGE'  
    SSFolder = 'Datafiles_SSUSI'
    GUFolder = 'Datafiles_GUVI'

    #Set filepaths to None initially.
    IMfp = None 
    SGfp = None
    
    #If there is a manual file input, set filepath for IMAGE and/or SSUSI/GUVI 
    if len(mfl) != 0:
        for fp in mfl: #fp = filepath
            if fp.split('.')[-1] == 'cdf':
                IMfp = fp
            elif fp.split('.')[-1] == 'nc':
                SGfp = fp
    
    #At this point there will be both IMfp and SGfp, one of either, or none.
    #Below codes will fill up whatever missing until there are filepaths for both.
    
    #If no IMAGE filepath yet, search for nearest datetime match.    
    if IMfp == None:
        IMfp = SP.IMAGE_filesearch(IMFolder, IMtype, dtime)
    
    #If no SSUSI/GUVI filepath yet, search for nearest datetime match.
    if SGfp == None:
        if SGind == 0:
            namestr = 'GUVI'
            GUtype = 'edr-aur'
            SGfp = SP.GUVI_filesearch(GUFolder, GUtype, dtime)
        elif SGind == 1:
            namestr = 'SSUSI'
            SStype = 'edr-aur'
            SSinstr = 'F16' #only one available during IMAGE lifetime
            SGfp = SP.SSUSI_filesearch(SSFolder, SStype, SSinstr, dtime)

    #Extract IMAGE data first.
    DataID = IMfp.split('/')[-1].split('__')[0]
    IMData = IMD.data(IMfp, DataID)
    
    #If there is only one datetime in IMData, then there is actually no data available because of filled value.
    if len(IMData.datetimes) == 1:
        print('Only one datetime exists in hourly file. Likely numerical filler barring incredible exception.')
        print(IMData.datetimes[0])
        return
    else:
        #Loop over datetimes inside IMAGE file
        if dtime not in IMData.datetimes:
            #If there is no match with footprint datetime, end function.
            for item in IMData.datetimes:
                if type(item) is not dt.datetime:
                    print('At least one errant datetime data in file which may cause no match')
                    return 
            
            #If there is match, transform dtime to nearest element in self.datetimes.
            diff_arr = abs(dtime-IMData.datetimes)
            min_pos = np.where(diff_arr == min(diff_arr))
            dtime = IMData.datetimes[min_pos][0]
            print('Producing closest match IMAGE figure')
        
        #Obtain IMAGE intensity data.
        im = IMData.image[IMData.datetimes==dtime,:,:]
        im = im[0,:,:]
        hemis = IMData.hemis[IMData.datetimes==dtime]
            
    #Plot IMAGE data in grayscale
    fig,ax = plt.subplots(1,1,figsize=(16,16),squeeze=False)
    
    if IMtype == 'WIC': #Set max-mim standard ranges
        cmin = 0
        cmax = 3500
        cmap_str = 'binary' #Dark where there is higher intensity
        title_type = 'WIC'
    elif IMtype == 'S12':
        cmin = 0
        cmax = 1200
        cmap_str = 'binary' 
        title_type = 'SI-12'
        
    #Set color mode.
    if mode == 'light':
        ctick = 'k'
        tchoice = False
    elif mode == 'dark':
        ctick = 'w'
        tchoice = True
    else:
        raise Exception('Incorrect Colour Mode')
    
    #Set x-y axes properties.
    xx, yy = np.meshgrid(IMData.XD1,IMData.YD2)
    ax[0,0].pcolormesh(xx,yy,im,vmin=cmin,vmax=cmax,cmap=cmap_str)
    ax[0,0].minorticks_on() #set minor axes locator
    ax[0,0].xaxis.set_minor_locator(AutoMinorLocator(5))
    ax[0,0].yaxis.set_minor_locator(AutoMinorLocator(5))
    ax[0,0].set_xlabel('i (km)',fontfamily='sans-serif',fontsize=32,labelpad=2)
    ax[0,0].set_ylabel('j (km)',fontfamily='sans-serif',fontsize=32,labelpad=0)
    ax[0,0].set_aspect('equal',adjustable='box')
    ax[0,0].tick_params(axis='both', which='major', direction='inout', width=1, length=15, labelsize=24)
    ax[0,0].tick_params(axis='both', which='minor', direction='inout', width=1, length=5, labelsize=12)
    ax[0,0].title.set_text(title_type + ' ' + dtime.strftime('%H:%M:%S'))
    ax[0,0].title.set_fontsize(48)
    ax[0,0].title.set_fontfamily('sans-serif')
    
    #Set axes color properties.
    ax[0,0].spines['bottom'].set_color(ctick) 
    ax[0,0].spines['top'].set_color(ctick)
    ax[0,0].spines['left'].set_color(ctick)
    ax[0,0].spines['right'].set_color(ctick)
    ax[0,0].tick_params(axis='both',which='both',colors=ctick) 
    ax[0,0].xaxis.label.set_color(ctick) 
    ax[0,0].yaxis.label.set_color(ctick)
    ax[0,0].title.set_color(ctick)
    
    
    #Obtain SG data for corresponding hemisphere to IMAGE
    SGData = SGD.SGExtract(SGfp)
    EDRdict = SGData.EDR_Aurora_Data
    MagLat = EDRdict['LATITUDE_GEOMAGNETIC_GRID_MAP'][:,:]
    MLT = EDRdict['MLT_GRID_MAP'][:,:]
    
    #Threshold each channel to remove (some) background, then sum up channel radiances.
    Int = np.zeros_like(MagLat)
    for ch in SGChan:
        if hemis == 0: #Select data for northern or southern polar cap.
            Int_ch = EDRdict['DISK_RADIANCEDATA_INTENSITY_NORTH'][ch-1,:,:]
        elif hemis == 1:
            Int_ch = EDRdict['DISK_RADIANCEDATA_INTENSITY_SOUTH'][ch-1,:,:]
        
        thres_min  = 100 #Set minimum threshold for plotting, specifically case of combined channels
        
        #Radiance threshold for each scan channel.
        if ch == 5: #not every channel is calibrated the same, so thresholds acan be different
            thres = 100
            if thres < thres_min: #If threshold is actually lower than minimum above, set a new minimum
                thres_min = thres 
        elif ch == 4:
            thres = 100
            if thres < thres_min:
                thres_min = thres
        elif ch == 3:
            thres = 100
            if thres < thres_min:
                thres_min = thres
        else:
            thres = thres_min
        
        #Sum up intensities across channels.
        Int_ch[Int_ch<=thres] = 0 
        Int += Int_ch
    
    #Convert SSUSI/GUVI AACGM latitudes and MLT into IMAGE Cartesian X-Y
    ang = MLT*np.pi/12
    ang = ang.astype(np.float64) #SSUSI/GUVI give data in float-64.
    xarr = 111*(90-abs(MagLat))*np.sin(ang) #scale*colat*lon
    yarr = -1*111*(90-abs(MagLat))*np.cos(ang)
    
    #Overlay SSUSI/GUVI data onto IMAGE.
    SGfig = ax[0,0].scatter(xarr,yarr,c=Int,s=1,norm=colors.LogNorm(vmin=thres_min,vmax=10000),
               cmap='turbo',zorder=1)
    
    #Set color bar and its color properties.
    cbar0 = fig.colorbar(SGfig,ax=ax[0,0],label='{} Radiance (R)'.format(namestr),pad=0.01)
    cbar0.outline.set_color(ctick)
    cbar0.ax.axes.tick_params(axis='y',which='both',colors=ctick,labelsize=32) 
    cbar0.ax.axes.yaxis.label.set_color(ctick)
    cbar0.ax.axes.yaxis.label.set_fontsize(34)
    
    #Adjust subplot margins.
    fig.subplots_adjust(left=0.07,right=0.97,bottom=0.15,top=0.85,wspace=0.0,hspace=0.0)
    
    #Save Figure
    if save == True:
        dtstr = dtime.strftime('%Y-%m-%d_%H-%M-%S')
        savename = 'SGonIM_{}.{}'.format(dtstr,saveformat)
        plt.savefig(savename, dpi=96, format=saveformat, transparent=tchoice)
        return fig,ax,hemis
    else:
        return fig,ax,hemis



"""
Example Script
"""
#Obtain SSUSI/GUVI on IMAGE figure.
dtime = dt.datetime(2002,3,18,14,55,0)
IMtype = 'S12'
SGind = 0 #0 = GUVI, 1 = SSUSI
SGChan = [3,4]
mfl = []
mode = 'light'
save = True
saveformat = 'png'
fig,ax,hemis = SGonIM(dtime,IMtype,SGind,SGChan,mfl,mode,save,saveformat) 


#Obtain Cluster footprint at datetime for further overlay.
dt_ft = [dtime] #Make datetime a list to use in footprint class.
sclist = ['C1']
FTFolder = 'Cluster_FT_T96_Maunder'
FT = FTOverlay(dt_ft,sclist,FTFolder)

#Set up correct hemisphere for SSUSI/GUVI, then set overlay figure to save.
hno_arr = np.array([hemis])
dt_union = np.array(dt_ft)
final_save = True
savefolder = 'IMSG_T96'
data_record = False

#Plot AACGM footprint coordinates.
fig,ax = OP.IMTSOverlay(fig,ax,IMtype,hno_arr,dt_union,FT.dt_dict,FT.agm_latN_dict,
                                 FT.agm_mltN_dict,FT.rhoN_dict,FT.agm_latS_dict,FT.agm_mltS_dict,
                                 FT.rhoS_dict,final_save,savefolder,mode,saveformat,data_record)



    