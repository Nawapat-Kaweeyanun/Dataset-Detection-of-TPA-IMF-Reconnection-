"""
Edition Date: 2025-August-29
@author: Nawapat Kaweeyanun
"""

"""
Objective: Plot IMAGE auroral images as a timeseries over an input set of datetimes

Prerequisite: A folder containing IMAGE CDF data files spanning the datetime range
"""
import os
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import SearchPlotter as SP
import IMAGE_Data as IMD
import csv


def timeseries_plot(dt_list,ncol,data_type,savefolder=None,save=False,mode='light',saveformat='png',data_record=False):
   
    #From a list of datetimes, determine number of subplots
    no_subplots = len(dt_list)
    
    #Determine number of rows needed from number of column given
    nrow = int(np.ceil(no_subplots/ncol))
    
    #Set up hemisphere list with same dimension as datetime list.
    hno_arr = []
    
        
    #Make figure
    figsize_ratio = nrow/ncol
    fig,ax = plt.subplots(nrow,ncol,figsize=(16,16*figsize_ratio),squeeze=False)
        
    if data_type == 'WIC': #Set intensity range for WIC and SI-12 (or S12)
        cmin = 1000
        cmax = 2000
        cmap_str = 'seismic'
        title_type = 'WIC'
    elif data_type == 'S12':
        cmin = 0
        cmax = 1200
        cmap_str = 'seismic'
        title_type = 'SI-12'
    
    #Set color modes
    if mode == 'light':
        ctick = 'k'
        tchoice = False
    elif mode == 'dark':
        ctick = 'w'
        tchoice = True
    else:
        raise Exception('Incorrect Colour Mode')
        
    #Plot one subplot at a time    
    for i in np.arange(no_subplots):
            
        #Obtain datetime
        dtime = dt_list[i]
            
        #Utilise search function to identify correct cdf file for this datetime
        IMFolder = 'Datafiles_IMAGE'
        filepath = SP.IMAGE_filesearch(IMFolder,data_type,dtime)
        
        #Extract data from file
        DataID = filepath.split('/')[-1].split('__')[0]
        IMData = IMD.data(filepath, DataID)
        
        #From all datetimes inside file, match closest datetime to input and obtain data
        #Account for cases when only one datetime exists in IMAGE file.
        if len(IMData.datetimes) == 1:
            
            if dtime not in IMData.datetimes:
                print('Datetime input not match one(s) in file')
                if type(IMData.datetimes[0]) is not dt.datetime:
                    print('No match because errant datetime data in file')
                    return
                else:
                    print('Producing closest match IMAGE figure')
                    dtime = IMData.datetimes[0] 
                    
                
            im = IMData.image[:,:]
            if IMData.hemis == 0:
                hno = 0
            elif IMData.hemis == 1:
                hno = 1
            else:
                hno = np.nan
            
        else:
    
            if dtime not in IMData.datetimes:
                print('Datetime input not match one(s) in file')
                for item in IMData.datetimes:
                    #Some files contains non-datetime item their datetime list ('corrupted')
                    if type(item) is not dt.datetime: 
                        print('At least one errant datetime data in file which may cause no match')
                        return 
                    
                #If there is no corrupted datetime, transform dtime to nearest element in self.datetimes.
                print('Producing closest match IMAGE figure')
                diff_arr = abs(dtime-IMData.datetimes)
                min_pos = np.where(diff_arr == min(diff_arr))
                dtime = IMData.datetimes[min_pos][0]
        
            im = IMData.image[IMData.datetimes==dtime,:,:]
            im = im[0,:,:]
            hemis = IMData.hemis[IMData.datetimes==dtime]
            if hemis == 0:
                hno = 0
            elif hemis == 1:
                hno = 1
            else:
                hno = np.nan
                    
        #Find subplot coordinate for this datetime. 
        #Plot from left to right, top to bottom
        row_ind = int(np.floor(i/ncol))
        col_ind = np.remainder(i,ncol)
        
        #Append hemisphere list
        hno_arr.append(hno)
        
        #Set axis and title fontsizes depending on subplot dimension
        if nrow != 1 and ncol != 1:
            axis_fontsize = 24
            title_fontsize = 24
        elif nrow != 1 and ncol == 1:
            axis_fontsize = 48
            title_fontsize = 48
        elif nrow == 1 and ncol != 1:
            axis_fontsize = 14
            title_fontsize = 16
        elif nrow == 1 and ncol == 1:
            axis_fontsize = 32
            title_fontsize = 48
        
        
        #Create X-Y meshgrid
        xx, yy = np.meshgrid(IMData.XD1,IMData.YD2)
        
        #Plot IMAGE subplot
        m = ax[row_ind,col_ind].pcolormesh(xx,yy,im,vmin=cmin,vmax=cmax,cmap=cmap_str)
        
        #Set axes properties
        ax[row_ind,col_ind].minorticks_on() 
        ax[row_ind,col_ind].xaxis.set_minor_locator(AutoMinorLocator(5))
        ax[row_ind,col_ind].yaxis.set_minor_locator(AutoMinorLocator(5))
        ax[row_ind,col_ind].set_xlabel('i (km)',fontsize=axis_fontsize,fontfamily='sans-serif',labelpad=2)
        ax[row_ind,col_ind].set_ylabel('j (km)',fontsize=axis_fontsize,fontfamily='sans-serif',labelpad=0) 
        ax[row_ind,col_ind].set_aspect('equal',adjustable='box')
        ax[row_ind,col_ind].tick_params(axis='both', which='major', direction='inout', width=1, length=15, labelsize=18)
        ax[row_ind,col_ind].tick_params(axis='both', which='minor', direction='inout', width=1, length=5, labelsize=18)
        
        #Add time title for each subplot
        ax[row_ind,col_ind].title.set_text(title_type + ' ' + dtime.strftime('%H:%M:%S'))
        ax[row_ind,col_ind].title.set_fontsize(title_fontsize)
        ax[row_ind,col_ind].title.set_fontfamily('sans-serif')
        
        
        #Set colors
        ax[row_ind,col_ind].spines['bottom'].set_color(ctick) #Set axes color
        ax[row_ind,col_ind].spines['top'].set_color(ctick)
        ax[row_ind,col_ind].spines['left'].set_color(ctick)
        ax[row_ind,col_ind].spines['right'].set_color(ctick)
        ax[row_ind,col_ind].tick_params(axis='both',which='both',colors=ctick) #Set tick color
        ax[row_ind,col_ind].xaxis.label.set_color(ctick) #Set x-axis ticklabel color
        ax[row_ind,col_ind].yaxis.label.set_color(ctick) #Set y-axis ticklabel color
        ax[row_ind,col_ind].title.set_color(ctick)
        
        #Write csv file for each subplot
        if data_record==True:
            dtstr = dtime.strftime('%Y-%m-%dT%H-%M-%SZ')
            csv_filename = '_'.join(['IMAGE',data_type,dtstr]) + '.csv'
            with open(csv_filename,'w', newline = '') as csvfile:
                csvwriter = csv.writer(csvfile)        
                for i in np.arange(len(IMData.XD1)):
                    row = [IMData.XD1[i],IMData.YD2[i]]
                    row = np.concatenate((row,im[i,:]))
                    csvwriter.writerow(row)
                csvfile.close()

        
    #Remove all subplots not used
    for j in np.arange(no_subplots,ncol*nrow): #All indices between no. plots filled and no. created.
        row_ind = int(np.floor(j/ncol))
        col_ind = np.remainder(j,ncol)
        ax[row_ind,col_ind].remove()
        
    #Set margins and add color bar axis, whose size depends on subplot dimension
    if no_subplots > 1:
        if nrow != 1 and ncol != 1:
            fig.subplots_adjust(left=0.09,right=0.88,bottom=0.01,top=0.99,wspace=0.3,hspace=0.0) #For 2X2
            cbar_ax = fig.add_axes([0.90,0.095,0.02,0.81]) #[hori_start,verti_start,width,height] as fraction of image, for 2x2
            cb = plt.colorbar(m,cax=cbar_ax)
            cb.ax.tick_params(labelsize=24)
            cb.set_label(label='Intensity (R)',size=24,family='sans-serif')
            cb.ax.xaxis.set_tick_params(pad=-60)
        elif nrow != 1 and ncol == 1: 
            fig.subplots_adjust(left=0.15,right=0.83,bottom=0.01,top=0.99,wspace=0.0,hspace=0.5) #For 1X4 
            cbar_ax = fig.add_axes([0.85,0.015,0.05,0.97]) 
            cb = plt.colorbar(m,cax=cbar_ax)
            cb.ax.tick_params(labelsize=32)
            cb.ax.xaxis.set_tick_params(pad=-60)
            cb.set_label(label='Intensity (R)',size=32,family='sans-serif')
        elif nrow == 1 and ncol != 1:
            fig.subplots_adjust(left=0.09,right=0.88,bottom=0.01,top=0.99,wspace=0.6,hspace=0.0) #For 1x4 Fig. 2A-2D and Fig. 2E-2H
            cbar_ax = fig.add_axes([0.89,0.095,0.02,0.81]) 
            cb = plt.colorbar(m,cax=cbar_ax)
            cb.ax.tick_params(labelsize=24)
            cb.ax.xaxis.set_tick_params(pad=-60)
            cb.set_label(label='Intensity (R)',size=24,family='sans-serif')
    elif no_subplots == 1:
        fig.subplots_adjust(left=0.03,right=0.95,bottom=0.01,top=0.99,wspace=0.10,hspace=0.0) #for 1X1
        cbar_ax = fig.add_axes([0.86,0.14,0.02,0.72])
        cb = plt.colorbar(m,cax=cbar_ax)
        cb.ax.tick_params(labelsize=32)
        cb.ax.xaxis.set_tick_params(pad=-60)
        cb.set_label(label='Intensity (R)',size=32,family='sans-serif')
        
    #Color bar colour settings
    cb.ax.axes.tick_params(axis='y',which='both',colors=ctick) #Set colorbar tick color
    cb.ax.axes.yaxis.label.set_color(ctick) #Set colorbar label color
    
    #Save file
    if save == True: #If wish to save
        dt_startstr = dt_list[0].strftime('%Y%m%d%H%M%S')
        dt_endstr = dt_list[-1].strftime('%Y%m%d%H%M%S')
        
        if savefolder == None:
            savefolder = os.getcwd() #If savefolder not given, use current
         
        name = 'IMAGE_Timeseries_{}_{}_{}.{}'.format(data_type,dt_startstr,dt_endstr,saveformat)
        savepath = '{}/{}/{}/{}/'.format(savefolder,data_type,dt_list[0].year,dt_list[0].month)
        os.makedirs(savepath,exist_ok=True)
        plt.savefig(savepath + name,dpi=300,format=saveformat,transparent=tchoice)
        plt.close()
        hno_arr = np.array(hno_arr) #for now
        return fig,ax,hno_arr
    elif save == False: #I doesn't wish to save
        plt.close()
        hno_arr = np.array(hno_arr)
        return fig,ax,hno_arr


"""
Example Script 
"""

"""
ncol = 4
data_type = 'WIC'
saveformat = 'svg'
savefolder = 'IMAGE_Timeseries'
save = True
dt_list = [dt.datetime(2002,3,18,14,15,0),
           dt.datetime(2002,3,18,14,30,0),
           dt.datetime(2002,3,18,14,45,0),
           dt.datetime(2002,3,18,15,0,0)]
#dt_list = [dt.datetime(2002,3,16,6,30,0)]
timeseries_plot(dt_list,ncol,data_type,savefolder,save,saveformat)
"""

