"""
Edition Date: 2025-September-02
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
from string import ascii_uppercase as ALC

    
def IMTSOverlay(IMfig,IMax,data_type,hno_arr,dt_union,dt_dict,latN_dict,mltN_dict,rhoN_dict,
                latS_dict,mltS_dict,rhoS_dict,save=False,savefolder=None,mode='light',saveformat='png',
                data_record=False):
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
    for i in np.arange(len(dt_union)):
        xdict.update({ALC[i]:{}})
        ydict.update({ALC[i]:{}})
    
        #Loop ver spacecraft in list.
        for key in dt_dict.keys():
            
            #If there is no data for this datetime, continue
            if np.size(np.where(dt_dict[key] == dt_union)[0][0]) == 0:
                continue
            
            #Choose northern or southern polar cap depending on data, then set footprint location
            if hno_arr[i] == 0: #If IMAGE over N cap
                if rhoN_dict[key][i] < 100+6371.2: #Footprint less than 100 km above surface
                    lat = latN_dict[key][i] #In deg
                    mlt = mltN_dict[key][i]
                    mark = 'o' #Default dot
                else: #Footprint never returned to ionosphere, flip to use the other side
                    lat = -1*latS_dict[key][i] #In deg
                    mlt = mltS_dict[key][i] 
                    mark = 'X' #Cross marker
            elif hno_arr[i] == 1: #If IMAGE over S cap
                if rhoS_dict[key][i] < 100+6371.2:
                    lat = latS_dict[key][i]
                    mlt = mltS_dict[key][i]
                    mark = 'o'
                else:
                    lat = -1*latN_dict[key][i]
                    mlt = mltN_dict[key][i]
                    mark = 'X'
                    
            else: #If error
                print('Cannot determine hemisphere for datetime {}. No footprint overlaid'.format(str(dt_union[i])))
                continue
            
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
            time_label = dt_union[i].strftime("%H:%M")
            
            #Add time label, if used.
            #IMax[row_ind,col_ind].text(x+label_offset_x,y+label_offset_y,time_label,
            #                           c=Colour_dict[key],fontsize=label_fontsize,zorder=2)
            
            #Add subplot letter label.
            IMax[row_ind,col_ind].text(text_xi,text_yi,'({})'.format(ALC[i]),fontsize=label_fontsize, fontfamily='sans-serif',
                                       horizontalalignment='center', verticalalignment='center',
                                       bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
            
            #Special 'I' subplot for IMAGE with SSUSI/GUVI overlay.
            #IMax[row_ind,col_ind].text(text_xi,text_yi,'(I)',fontsize=label_fontsize, fontfamily='sans-serif',
            #                           horizontalalignment='center', verticalalignment='center',
            #                           bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
        
    #Save figure        
    if save == True:
        figname = 'IMAGE_Timeseries_Footprint_{}_{}.{}'.format(data_type,dt_union[0].strftime('%Y-%m-%d_%H-%M-%S'),saveformat) 
        savedir = savefolder + '/{}/{}/{}/'.format(data_type,dt_union[0].year,dt_union[0].month)
        os.makedirs(savedir,exist_ok=True)
        plt.savefig(savedir + figname,dpi=200,format=saveformat,transparent=tchoice)
        plt.close()
        
        #Record footprint location to a csv file.
        if data_record == True:
            date_str = str(dt_union[0]).split(' ')[0]
            csv_filename = '_'.join(['ClusterFootprint',data_type,date_str]) + '.csv'
            with open(csv_filename,'w', newline = '') as csvfile:
                csvwriter = csv.writer(csvfile)
                        
                for i in np.arange(len(dt_union)):
                    row = [dt_union[i]]
                    for key in dt_dict.keys():
                        row.append(xdict[ALC[i]][key])
                        row.append(ydict[ALC[i]][key])
                    csvwriter.writerow(row)
            csvfile.close()
        
        return IMfig,IMax

    else:
        return IMfig,IMax
        

def SGOverlay(SGfig,SGax,ft_dt,ft_latN,ft_mltN,ft_rhoN,ft_latS,ft_mltS,
              ft_rhoS,satellite,SGinstr,Chan,OrbNo,switch=True,mode='light',
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
    
    #Plot latitude points.
    #Note: SSUSI data has southern cap latitude as positive values. The negative latitude axis label is simply labelling.
    #Convert S cap latitude to positive to be plotted.
    
    for key in ft_latN.keys():
        if ft_rhoN[key] < 100+6371.2: #False if nan
            SGax[0].scatter(np.pi+(ft_mltN[key]*np.pi/12),ft_latN[key],c=Colour_dict[key],marker='o',s=100)
            SGax[0].text(np.pi+((ft_mltN[key]+label_offset_mlt)*np.pi/12),ft_latN[key]+label_offset_lat,
                         time_label,c='magenta',fontsize=label_fontsize)
        elif switch == True:
            SGax[0].scatter(np.pi+(ft_mltS[key]*np.pi/12),abs(ft_latS[key]),c=Colour_dict[key],marker='X',s=100)
            SGax[0].text(np.pi+((ft_mltS[key]+label_offset_mlt)*np.pi/12),abs(ft_latS[key])+label_offset_lat,
                         time_label,c='magenta',fontsize=label_fontsize)
        else:
            continue
        
        
        if ft_rhoS[key] < 100+6371.2:    
            SGax[1].scatter(np.pi+(ft_mltS[key]*np.pi/12),abs(ft_latS[key]),c=Colour_dict[key],marker='o',s=100)
            SGax[1].text(np.pi+((ft_mltS[key]+label_offset_mlt)*np.pi/12),abs(ft_latS[key])+label_offset_lat,
                         time_label,c='magenta',fontsize=label_fontsize)
        elif switch == True:
            SGax[1].scatter(np.pi+(ft_mltN[key]*np.pi/12),ft_latN[key],c=Colour_dict[key],marker='X',s=100)
            SGax[1].text(np.pi+((ft_mltN[key]+label_offset_mlt)*np.pi/12),ft_latN[key]+label_offset_lat,
                         time_label,c='magenta',fontsize=label_fontsize)
        else:
            continue
    
    
    #Save figure
    if save == True:
        
        chan_str = '+'.join([str(ch) for ch in Chan])
        ft_timestr = dt.datetime.strftime(ft_dt, '%H%M%S')
        
        #Different file locations for SSUSI and GUVI.
        if satellite == 'SSUSI':
            figname = 'SSUSI_T96_{}_{}_{}_{}.png'.format(SGinstr,chan_str,OrbNo,ft_timestr) 
            savedir = savefolder + '/{}/{}/{:03d}/Channel {}/'.format(SGinstr,ft_dt.year,ft_dt.timetuple().tm_yday,chan_str)
        elif satellite == 'GUVI':
            figname = 'GUVI_T96_{}_{}_{}.png'.format(chan_str,OrbNo,ft_timestr)
            savedir = savefolder + '/{}/{:03d}/Channel {}/'.format(ft_dt.year,ft_dt.timetuple().tm_yday,chan_str)
        
        os.makedirs(savedir,exist_ok=True)
        plt.savefig(savedir + figname,dpi=100,format=saveformat,transparent=tchoice)
        plt.close()
        return SGfig,SGax
    
    else:
        return SGfig,SGax