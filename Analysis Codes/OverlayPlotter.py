"""
Edition Date: 2025-October-21
@author: Nawapat Kaweeyanun
"""

"""
Objective: Functions for plotting T96 footprints onto IMAGE or SSUSI/GUVI snapshots.
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import csv
import aacgmv2
from string import ascii_uppercase as ALC

    
def IMTSOverlay(IMfig,IMax,data_type,hno_arr,FTdt_union,FTdt_dict,latN_dict,mltN_dict,rhoN_dict,
                latS_dict,mltS_dict,rhoS_dict,IMdt_union,HiRes=False,save=False,savefolder=None,
                mode='light',saveformat='png',data_record=False):
    #Method for plotting footprint for IMAGE timeseries
    #Important: dt_union must be the same array used to create IMAGE timeseries!
    #Coordinates must correspond to datetimes
    
    #Set colour for each spacecraft
    Colour_dict = {}
    Colour_dict.update({'C1':'gold'})
    Colour_dict.update({'C2':'r'})
    Colour_dict.update({'C3':'magenta'})
    Colour_dict.update({'C4':'b'})
    
    #Set colour modes
    if mode == 'light':
        tchoice = False
    elif mode == 'dark':
        tchoice = True
    else:
        raise Exception('Incorrect Colour Mode')
    
    #Open IMAGE figure from input.
    plt.figure(IMfig)
    
    #Dictionaries for final x-y positions to plot
    xdict = {} 
    ydict = {}
 
    #Loop over footprint datetime array 
    for i in np.arange(len(FTdt_union)):
        xdict.update({ALC[i]:{}})
        ydict.update({ALC[i]:{}})
    
        #Loop ver spacecraft in list.
        for key in FTdt_dict.keys():
            
            #If there is no data for this datetime, continue
            if np.size(np.where(FTdt_dict[key] == FTdt_union)[0][0]) == 0:
                continue
            
            #Choose northern or southern polar cap depending on data, then set footprint location
            if hno_arr[i] == 0: #If IMAGE over N cap
                if rhoN_dict[key][i] < 200+6371.2: #Footprint less than 100 km above surface
                    lat = latN_dict[key][i] #In deg
                    mlt = mltN_dict[key][i]
                    mark = 'o' #Default dot
                else: #Footprint never returned to ionosphere, flip to use the other side
                    lat = -1*latS_dict[key][i] #In deg
                    mlt = mltS_dict[key][i] 
                    mark = 'X' #Cross marker
            elif hno_arr[i] == 1: #If IMAGE over S cap
                if rhoS_dict[key][i] < 200+6371.2:
                    lat = latS_dict[key][i]
                    mlt = mltS_dict[key][i]
                    mark = 'o'
                else:
                    lat = -1*latN_dict[key][i]
                    mlt = mltN_dict[key][i]
                    mark = 'X'
                    
            else: #If error
                print('Cannot determine hemisphere for datetime {}. No footprint overlaid'.format(str(FTdt_union[i])))
                continue
            
            #If IMAGE data is from ECLAT.
            if HiRes == False: 
                #Convert footprint from AACGM to Cartesian
                lon = mlt*np.pi/12 #In radians
                colat = (90-abs(lat)) #In deg
                x = 111*colat*np.sin(lon)
                y = 111*-1*colat*np.cos(lon)
                xdict[ALC[i]].update({key:x}) #Update footprint dictionaries.
                ydict[ALC[i]].update({key:y})
                
                #Set footprint marker size, time label position offset (if used), and subplot label position
                if np.size(IMax) == 1: #Enlarge sizes if there is only one subplot.
                    msize = 500
                    label_offset_x = 200
                    label_offset_y = 300
                    label_fontsize = 48
                    text_xi = -3900
                    text_yi = 4000
                else:
                    msize = 100
                    label_offset_x = 200
                    label_offset_y = 300
                    label_fontsize = 14
                    text_xi = -3600
                    text_yi = 3800
                    
                
                #Identify correct subplot for each footprint.
                [nrow,ncol] = np.shape(IMax)
                row_ind = int(np.floor(i/ncol))
                col_ind = np.remainder(i,ncol)
                
                #Overlay footprint on IMAGE subplot.
                IMax[row_ind,col_ind].scatter(x,y,c=Colour_dict[key],marker=mark,s=msize)
                time_label = FTdt_union[i].strftime("%H:%M")
                
                #Add time label, if used.
                IMax[row_ind,col_ind].text(x+label_offset_x,y+label_offset_y,time_label,
                                           c=Colour_dict[key],fontsize=label_fontsize,zorder=2)
                
                #Add subplot letter label.
                IMax[row_ind,col_ind].text(text_xi,text_yi,'({})'.format(ALC[i]),fontsize=label_fontsize, fontfamily='sans-serif',
                                           horizontalalignment='center', verticalalignment='center',
                                           bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                
                #Special 'I' subplot for IMAGE with SSUSI/GUVI overlay.
                #IMax[row_ind,col_ind].text(text_xi,text_yi,'(I)',fontsize=label_fontsize, fontfamily='sans-serif',
                #                           horizontalalignment='center', verticalalignment='center',
                #                           bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
            
            #If use high-resolution IMAGE data, no need to convert AACGM coordinate to ECLAT Cartesian
            elif HiRes == True:
                
                #Set dot and label properties
                time_label = FTdt_union[i].strftime('%H:%M')
                label_offset_mlt = -0.3
                label_offset_lat = -4
                msize = 40 #75 #Fig2A-L #40 Fig2M
                label_fontsize = 12
                
                [nrow,ncol] = np.shape(IMax)
                row_ind = int(np.floor(i/ncol))
                col_ind = np.remainder(i,ncol)
                IMax[row_ind,col_ind].scatter(np.pi+(mlt*np.pi/12),lat,c=Colour_dict[key],marker='o',s=msize,zorder=2)
                IMax[row_ind,col_ind].text(np.pi+((mlt+label_offset_mlt)*np.pi/12),lat+label_offset_lat,time_label,c=Colour_dict[key],fontsize=label_fontsize,zorder=2)
                
                #Update footprint dictionary with mlt in 'x' and latitude in 'y' components.
                xdict[ALC[i]].update({key:mlt})
                ydict[ALC[i]].update({key:lat})
        
    #Save figure        
    if save == True:
        figname = 'IMAGE_{}_{}_{}_Footprint_{}_{}.{}'.format(data_type,IMdt_union[0].strftime('%Y%m%d%H%M%S'),
                                                          IMdt_union[-1].strftime('%Y%m%d%H%M%S'),
                                                          FTdt_union[0].strftime('%Y%m%d%H%M%S'),
                                                          FTdt_union[-1].strftime('%Y%m%d%H%M%S'),saveformat) 
        savedir = savefolder + '/{}/{}/'.format(FTdt_union[0].year,FTdt_union[0].month)
        os.makedirs(savedir,exist_ok=True)
        plt.savefig(savedir + figname,dpi=200,format=saveformat,transparent=tchoice)
        plt.close()
        
        #Record footprint location to a csv file.
        if data_record == True:
            date_str = str(FTdt_union[0]).split(' ')[0]
            csv_filename = '_'.join(['ClusterFootprint',data_type,date_str]) + '.csv'
            with open(csv_filename,'w', newline = '') as csvfile:
                csvwriter = csv.writer(csvfile)
                        
                for i in np.arange(len(FTdt_union)):
                    row = [FTdt_union[i]]
                    for key in FTdt_dict.keys():
                        row.append(xdict[ALC[i]][key])
                        row.append(ydict[ALC[i]][key])
                    csvwriter.writerow(row)
            csvfile.close()
        
        return IMfig,IMax

    else:
        return IMfig,IMax
        

def SGOverlay(SGfig,SGax,ft_dt,ft_latN,ft_mltN,ft_rhoN,ft_latS,ft_mltS,
              ft_rhoS,sc,Chan,OrbNo,cap,switch=True,mode='light',
              save=False,savefolder=None,saveformat='png'):
    #Method for overlay Cluster footprint on SSUSI/GUVI auroral scans.
    #Both northern and southern footprints because SSUSI/GUVI plots contain both
    #Footprint locations are dictionaries with spacecraft as keys
    #ACCEPT ONLY ONE DATETIME. Call plot multiple times for multiple datetimes
    
    #Set colour for each spacecraft
    Colour_dict = {}
    Colour_dict.update({'C1':'gold'})
    Colour_dict.update({'C2':'r'})
    Colour_dict.update({'C3':'magenta'})
    Colour_dict.update({'C4':'b'})
    
    #Set Color mode
    if mode == 'light':
        tchoice = False
    elif mode == 'dark':
        tchoice = True
    else:
        raise Exception('Incorrect Colour Mode')
    
    #Open SSUSI/GUVI scan.
    plt.figure(SGfig)
    time_label = ft_dt.strftime('%H:%M')
    label_offset_mlt = -0.3
    label_offset_lat = -4
    label_fontsize = 16
    
    #Set plot parameters depending on whether plotting one (north or south cap only) or two subplots (both).
    if cap in ['N','S']:
        ax_ind = [0] #One subplot index
        dot_size = 25
        label_fontsize = 10
    elif cap == 'B':
        ax_ind = [0,1] #Two subplot indices
        dot_size = 100
        label_fontsize = 16
    
    #Plot footprint.
    for key in ft_latN.keys(): #Loop over Cluster list
       for ind in ax_ind:
           if ind == 0:
               #If northern cap, check for northern footprint first, if not then switch for conjugate footprint.
               #For both-caps plot, northern cap will always come first, so can condense some if-else branches.
               if cap in ['N','B']: 
                   if ft_rhoN[key] < 200+6371.2:
                       SGax[0,ind].scatter(np.pi+(ft_mltN[key]*np.pi/12),ft_latN[key],c=Colour_dict[key],marker='o',s=dot_size)
                       SGax[0,ind].text(np.pi+((ft_mltN[key]+label_offset_mlt)*np.pi/12),ft_latN[key]+label_offset_lat,
                                    time_label,c=Colour_dict[key],fontsize=label_fontsize)
                   elif switch == True:
                       SGax[0,ind].scatter(np.pi+(ft_mltS[key]*np.pi/12),abs(ft_latS[key]),c=Colour_dict[key],marker='X',s=dot_size)
                       SGax[0,ind].text(np.pi+((ft_mltS[key]+label_offset_mlt)*np.pi/12),abs(ft_latS[key])+label_offset_lat,
                                    time_label,c=Colour_dict[key],fontsize=label_fontsize)
               #If southern cap, check for southern footprint first, if not then switch for conjugate footprint.
               elif cap == 'S': 
                   if ft_rhoS[key] < 200+6371.2:
                       SGax[0,ind].scatter(np.pi+(ft_mltS[key]*np.pi/12),abs(ft_latS[key]),c=Colour_dict[key],marker='o',s=dot_size)
                       SGax[0,ind].text(np.pi+((ft_mltS[key]+label_offset_mlt)*np.pi/12),abs(ft_latS[key])+label_offset_lat,
                                    time_label,c=Colour_dict[key],fontsize=label_fontsize)
                   elif switch == True:
                       SGax[0,ind].scatter(np.pi+(ft_mltN[key]*np.pi/12),ft_latN[key],c=Colour_dict[key],marker='X',s=dot_size)
                       SGax[0,ind].text(np.pi+((ft_mltN[key]+label_offset_mlt)*np.pi/12),ft_latN[key]+label_offset_lat,
                                    time_label,c=Colour_dict[key],fontsize=label_fontsize)
           elif ind == 1:
               if ft_rhoS[key] < 200+6371.2:    
                   SGax[0,ind].scatter(np.pi+(ft_mltS[key]*np.pi/12),abs(ft_latS[key]),c=Colour_dict[key],marker='o',s=dot_size)
                   SGax[0,ind].text(np.pi+((ft_mltS[key]+label_offset_mlt)*np.pi/12),abs(ft_latS[key])+label_offset_lat,
                                    time_label,c=Colour_dict[key],fontsize=label_fontsize)
               elif switch == True:
                   SGax[0,ind].scatter(np.pi+(ft_mltN[key]*np.pi/12),ft_latN[key],c=Colour_dict[key],marker='X',s=dot_size)
                   SGax[0,ind].text(np.pi+((ft_mltN[key]+label_offset_mlt)*np.pi/12),ft_latN[key]+label_offset_lat,
                                    time_label,c=Colour_dict[key],fontsize=label_fontsize)
    
    
    #Save figure
    if save == True:
        
        chan_str = '+'.join([str(ch) for ch in Chan])
        ft_timestr = dt.datetime.strftime(ft_dt, '%H%M%S')
        
        #Different file locations for SSUSI and GUVI.
        if sc == 'TIMED':
            figname = 'GUVI_T96_{}_{}_{}_{}.png'.format(chan_str,OrbNo,ft_timestr,cap)
            savedir = savefolder + '/{}/{:03d}/Channel {}/'.format(ft_dt.year,ft_dt.timetuple().tm_yday,chan_str)
        else:
            figname = 'SSUSI_T96_{}_{}_{}_{}_{}.png'.format(sc,chan_str,OrbNo,ft_timestr,cap) 
            savedir = savefolder + '/{}/{}/{:03d}/Channel {}/'.format(sc,ft_dt.year,ft_dt.timetuple().tm_yday,chan_str)
        
        os.makedirs(savedir,exist_ok=True)
        plt.savefig(savedir + figname,dpi=100,format=saveformat,transparent=tchoice)
        plt.close()
        return SGfig,SGax
    
    else:
        return SGfig,SGax
    
def SCPath_Manual(filepath,dti,dtf,plotmode,HiRes,fig,ax,figmode,saveformat,save):
    #Plotter for DMSP/TIMED path from a SSC Locator file.
    
    #Open file and convert data into an array.
    with open(filepath) as csvfile:
        Data = csv.reader(csvfile,delimiter = ',')
        data_arr = []
        for row in Data:
            if len(data_arr) == 0: #First row.
                data_arr = np.array(row)
            else: #Then stack next rows under the first row.
                data_arr = np.vstack((data_arr,np.array(row)))
    
    #Use dictionary to store data with parameters as keys.
    DataDict = {}
    paralist = data_arr[0,:] #First row of file is header where parameter list resides.
    
    #Clean up some random characters that can appear.
    if 'Datetime' in paralist[0] and paralist[0] != 'Datetime': 
        paralist[0] = 'Datetime'
    
    #Update data dictionary. Turn datetime strings into datetimes and numeric strings into floats.
    for i in np.arange(len(paralist)):
        para_data = data_arr[1:,i]
        if 'Datetime' in paralist[i]: 
            para_data = np.array([dt.datetime.strptime(dtime,'%Y-%m-%dT%H:%M:%SZ') for dtime in para_data])
        else:
            para_data = np.array([float(item) for item in para_data])
        DataDict.update({paralist[i]:para_data})
    
    #Obtain spacecraft trajectory data from the data dictionary.
    dt_arr = DataDict['Datetime']
    geolat_arr = DataDict['Geolat']
    geolon_arr = DataDict['Geolon']
    
    #Filter out datetimes outside defined plotting windows
    #If starting datetime is before 1st datetime available in file (or final datetime after last file datetime),
    #then just get everything.
    low_lim = np.where(dt_arr >= dti)[0][0]
    high_lim = np.where(dt_arr > dtf)[0][0] #Not >= because of inclusive:exclusive nature of Python indices.
    dt_arr = dt_arr[low_lim:high_lim]
    geolat_arr = geolat_arr[low_lim:high_lim]
    geolon_arr = geolon_arr[low_lim:high_lim]

    #Convert path coordinate points from geocentric to aacgm.    
    agmlat_arr = []
    agmmlt_arr = []
    alt = 0
    for i in np.arange(len(dt_arr)):
        agm = aacgmv2.get_aacgm_coord(geolat_arr[i],geolon_arr[i],
                                      alt,dt_arr[i],method='GEOCENTRIC')
        agmlat_arr.append(agm[0])
        agmmlt_arr.append(agm[2])

    agmlat_arr = np.array(agmlat_arr)
    agmmlt_arr = np.array(agmmlt_arr)
    
    #Set path color properties depending on light or dark modes.
    if figmode.lower() == 'light':
        ctick = 'red'
        tchoice = False
    elif figmode.lower() == 'dark':
        ctick = 'red'
        tchoice = True
    else:
        raise Exception('Incorrect Colour Mode')
    
    #Separate plotmode, depend if applying path to IMAGE or SSUSI/GUVI figures.
    if plotmode == 'IM' and HiRes == False: #Original ECLAT IMAGE (one subplot only)
        #Convert path AACGM to ECLAT Cartesian
        lon = agmmlt_arr*np.pi/12
        colat = (90-abs(agmlat_arr))
        x_arr = 111*colat*np.sin(lon)
        y_arr = 111*-1*colat*np.cos(lon)

        #Annotate markers every 2 minutes.
        interval = 2
        x_mark = x_arr[0::2]
        y_mark = y_arr[0::2]
        dt_mark = dt_arr[0::2]
        
        #Plot path onto one IMAGE figure.
        plt.figure(fig)
        ax[0,0].plot(x_mark,y_mark,ctick,linewidth=2,marker='|',markevery=interval,zorder=1)
        for d in np.arange(len(dt_mark)):
            dtime_str = dt_mark[d].strftime('%M')
            ax[0,0].annotate(dtime_str,(x_mark[d],y_mark[d]),zorder=1,fontsize=24,color=ctick)
    
    elif plotmode == 'IM' and HiRes == True: #High-resolution IMAGE.
        #Annotate markers every 2 minutes.
        interval = 2
        agmmlt_mark = agmmlt_arr[0::2]
        agmlat_mark = agmlat_arr[0::2]
        dt_mark = dt_arr[0::2]
        
        #Plot AACGM path without Cartesian conversion.
        plt.figure(fig)
        ax[0,0].plot(np.pi+(agmmlt_mark*np.pi/12),agmlat_mark,c=ctick,ls='--',lw=1,marker='|',markevery=interval,zorder=3)
        for d in np.arange(len(dt_mark)):
            dtime_str = dt_mark[d].strftime('%M')
            ax[0,0].annotate(dtime_str,(np.pi+(agmmlt_mark[d]*np.pi/12),agmlat_mark[d]),zorder=4,fontsize=6,color=ctick)
    
    elif plotmode == 'SG': #SSUSI/GUVI. This has not been required so far. To be added later on.
        pass
    
    #Save figure.
    if save == True:
        date_str = dt_arr[0].strftime('%Y-%m-%d-%H-%M')
        figname = 'IMSG_Path_{}.{}'.format(date_str, saveformat)  #maybe change this to distinguish from CDAWeb Output
        plt.savefig(figname,dpi=200,format='png',transparent=tchoice)
        plt.close()
        return fig,ax
    else:
        return fig,ax
