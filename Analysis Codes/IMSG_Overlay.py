"""
Edition Date: 2025-October-29
@author: Nawapat Kaweeyanun
"""


"""
Objective: Overlay SSUSI/GUVI on top of IMAGE (greyscaled)

Rerequisite: Run from same directory as IMAGE and SSUSI/GUVI data directories.

Put all IDL files listed below in same directory
idlfilelist = ['wic20020771430.idl','wic20020771444.idl','wic20020771454.idl','wic20020771509.idl',
               's1320020771430.idl','s1320020771444.idl','s1320020771454.idl','s1320020771509.idl',
               's1220020771430.idl','s1220020771444.idl','s1220020771454.idl','s1220020771509.idl']

"""

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
from matplotlib.ticker import AutoMinorLocator
import SearchPlotter as SP
import IMAGE_Data as IMD
import SSUSI_GUVI_Data as SGD
import OverlayPlotter as OP
from Footprint_Aurora_Overlay import FTOverlay


def SGonIM(dtime,IMtype,IMHiRes=False,SGind=0,SGChan=[5],mfl=[],mode='light',cbar_on=True,save=False,saveformat='png'): #mfl = manual file input list
    #Function to overlay SSUSI/GUVI scan on grayscaled IMAGE auroral picture.
    #ACCEPT ONLY ONE DATETIME!

    if IMHiRes == False:
    
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
        if cbar_on == True:
            cbar0 = fig.colorbar(SGfig,ax=ax[0,0],label='{} Radiance (R)'.format(namestr),pad=0.01)
            cbar0.outline.set_color(ctick) #set colorbar outline color.
            cbar0.ax.axes.tick_params(axis='y',which='both',colors=ctick,labelsize=32) #set colorbar tick color
            cbar0.ax.axes.yaxis.label.set_color(ctick) #set colorbar label color
            cbar0.ax.axes.yaxis.label.set_fontsize(34)
        elif cbar_on == False:
            plt.tight_layout()
        
        #Adjust subplot margins.
        fig.subplots_adjust(left=0.07,right=0.97,bottom=0.15,top=0.85,wspace=0.0,hspace=0.0)
    
    elif IMHiRes == True: #HighRes IMAGE
        #Same principles as above but a new IMAGE file search method.
        
        SSFolder = 'Datafiles_SSUSI'
        GUFolder = 'Datafiles_GUVI'
        SGfp = None
        
        #Update SSUSI/GUVI file path if manual input exists.
        if len(mfl) != 0:
            for fp in mfl:
                if fp.split('.')[-1] == 'cdf': #IMAGE untouched.
                    pass
                elif fp.split('.')[-1] == 'nc':
                    SGfp = fp
        
        #Search for SSUSI/GUVI files if no manual inputs.
        if SGfp == None:
            if SGind == 0:
                namestr = 'GUVI'
                GUtype = 'edr-aur'
                SGfp = SP.GUVI_filesearch(GUFolder, GUtype, dtime)
            elif SGind == 1:
                namestr = 'SSUSI'
                SStype = 'edr-aur'
                SSinstr = 'F16'
                SGfp = SP.SSUSI_filesearch(SSFolder, SStype, SSinstr, dtime)
        
        idlfilelist = ['wic20020771430.idl','wic20020771444.idl','wic20020771454.idl','wic20020771509.idl',
                       's1320020771430.idl','s1320020771444.idl','s1320020771454.idl','s1320020771509.idl',
                       's1220020771430.idl','s1220020771444.idl','s1220020771454.idl','s1220020771509.idl']
        
        #Set IMAGE-related properties (color map limits)
        if IMtype == 'WIC':
            vmin = 500
            vmax = 5000
            title_name = 'WIC'
        elif IMtype == 'S12':
            vmin = 100
            vmax = 5000
            title_name = 'SI-12'
        elif IMtype == 'S13':
            vmin = 100
            vmax = 5000
            title_name = 'SI-13'
            
        #Obtain IMAGE data
        IHS = IMD.ImHighRes(idlfilelist,IMtype,[dtime])
        dtp_list = list(IHS.DataDict.keys())
        imp = IHS.DataDict[dtp_list[0]]['image'].ravel()
        mlatp = IHS.DataDict[dtp_list[0]]['mlat'].ravel()
        mlt = IHS.DataDict[dtp_list[0]]['mlt'].ravel()
        
        
        #Create IMAGE plot in grayscale.
        fig,ax = plt.subplots(1,1,subplot_kw={'projection': 'polar'},figsize=(4,4),squeeze=False)
        ax[0,0].scatter(np.pi+(mlt*np.pi/12),mlatp,c=imp,s=2,norm=colors.LogNorm(vmin=vmin,vmax=vmax),
                        cmap='binary',zorder=1)
        
        #Set angular axis properties
        ax[0,0].set_theta_zero_location('N')
        ax[0,0].set_rlim(90,40,5)
        xticks_loc = ax[0,0].get_xticks().tolist()
        ax[0,0].xaxis.set_major_locator(mticker.FixedLocator(xticks_loc))
        ax[0,0].set_xticklabels(['12','15','18','21','0','3','6','9'],fontsize=12)
        ax[0,0].spines['polar'].set_color(ctick) #set angular axis color
        
        #Set radial axis properties and subplot title
        yticks_loc = ax[0,0].get_yticks().tolist()
        if dtime.year in [2000,2001,2002,2003]:
            hemis = 0
            ax[0,0].yaxis.set_major_locator(mticker.FixedLocator(yticks_loc))
            ax[0,0].set_yticklabels(['','50','60','70','',''],fontsize=10)
            ax[0,0].set_title('IMAGE {} North\n{} (Grayscale)'.format(title_name,dtime.strftime('%H:%M')),color=ctick,fontsize=12,fontfamily='sans-serif')
        elif dtime.year in [2004,2005]:
            hemis = 1
            ax[0,0].yaxis.set_major_locator(mticker.FixedLocator(yticks_loc))
            ax[0,0].set_yticklabels(['','-50','-60','-70','',''],fontsize=10)
            ax[0,0].set_title('IMAGE {} South\n{} (Grayscale)'.format(title_name,dtime.strftime('%H:%M')),color=ctick,fontsize=12,fontfamily='sans-serif')
        ax[0,0].yaxis.label.set_color(ctick) #set y-axis ticklabel color
        ax[0,0].tick_params(axis='both',which='both',colors=ctick) #set tick color both axes
        ax[0,0].grid(axis='both',which='both',color=ctick) #set grid color
        
        #Obtain SG data in the same hemisphere as IMAGE.
        SGData = SGD.SGExtract(SGfp)
        EDRdict = SGData.EDR_Aurora_Data
        MagLat = EDRdict['LATITUDE_GEOMAGNETIC_GRID_MAP'][:,:]
        MLT = EDRdict['MLT_GRID_MAP'][:,:]
        
        #Sum up SSUSI/GUVI data cross channels as required
        Int = np.zeros_like(MagLat)
        for ch in SGChan:
            if hemis == 0:    
                Int_ch = EDRdict['DISK_RADIANCEDATA_INTENSITY_NORTH'][ch-1,:,:]
            elif hemis == 1:
                Int_ch = EDRdict['DISK_RADIANCEDATA_INTENSITY_SOUTH'][ch-1,:,:]
            thres_min  = 100 #final minimum threshold for plotting, specifically case of combined channels
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
            Int_ch[Int_ch<=thres] = 0 
            Int += Int_ch
    
        #Overlay SSUSI/GUVI plot onto the IMAGE plot.
        SGfig = ax[0,0].scatter(np.pi+(MLT*np.pi/12),MagLat,c=Int,s=1,norm=colors.LogNorm(vmin=thres_min,vmax=10000),
                   cmap='turbo',zorder=2)
        
        #Set color bar and figure layouts.
        if cbar_on == True:
            fig.subplots_adjust(left=0.10,right=0.72,bottom=0.05,top=0.83,wspace=0.0,hspace=0.0)
            cbar_ax = fig.add_axes([0.82,0.05,0.03,0.75]) #set axis location
            cbar0 = fig.colorbar(SGfig,cax=cbar_ax,label='{} Radiance (R)'.format(namestr),pad=0.01)
            cbar0.outline.set_color(ctick) #set colorbar outline color.
            cbar0.ax.axes.tick_params(axis='y',which='both',colors=ctick,labelsize=10) #set colorbar tick color
            cbar0.ax.axes.yaxis.label.set_color(ctick) #set colorbar label color
            cbar0.ax.axes.yaxis.label.set_fontsize(12)
        else:
            plt.tight_layout()
    
    
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
dtime = dt.datetime(2002,3,18,14,54,0)
IMtype = 'S12'
IMHiRes = True
SGind = 0 #0 = GUVI, 1 = SSUSI
SGChan = [3,4]
mfl = []
mode = 'dark'
cbar_on = False
save = False
saveformat = 'png'
fig,ax,hemis = SGonIM(dtime,IMtype,IMHiRes,SGind,SGChan,mfl,mode,cbar_on,save,saveformat)

#Obtain Cluster footprint at datetime for further overlay.
dt_ft = [dtime] #Make datetime a list to use in footprint class.
sclist = ['C1']
FTFolder = 'Cluster_FT_T96_Updated_127km' #Folder where Cluster footprint data is kept
FT = FTOverlay(dt_ft,sclist,FTFolder)

#Set up correct hemisphere for SSUSI/GUVI, then set overlay figure to save.
hno_arr = np.array([hemis])
FTdt_union = np.array(dt_ft)
IMdt_union = FTdt_union
save2 = False
savefolder = 'IMSG_T96'
data_record = False

#Plot AACGM footprint coordinates.
fig,ax = OP.IMTSOverlay(fig,ax,IMtype,hno_arr,FTdt_union,FT.dt_dict,FT.agm_latN_dict,
                                 FT.agm_mltN_dict,FT.rhoN_dict,FT.agm_latS_dict,FT.agm_mltS_dict,
                                 FT.rhoS_dict,IMdt_union,IMHiRes,save2,savefolder,mode,saveformat,data_record)

#Overlay TIMED/DMSP path.
SCfilepath = 'TIMED_Geo_18Mar02.csv' #target file containing TIMED location.
dti = dt.datetime(2002,3,18,14,49,0)
dtf = dt.datetime(2002,3,18,15,9,0)
final_save = True
fig,ax = OP.SCPath_Manual(SCfilepath,dti,dtf,'IM',IMHiRes,fig,ax,mode,saveformat,final_save)
    