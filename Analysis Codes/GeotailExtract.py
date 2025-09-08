
"""
Edition Date: 2025-August-29
@author: Nawapat Kaweeyanun
"""
"""
Objective: Extract & Plot Geotail Data from a downloaded file.

Prerequisite: LEP data (Editor A) downloaded from https://darts.isas.jaxa.jp/stp/geotail/ as a .txt.gz file in the same directory.
"""


import gzip
import numpy as np
import datetime as dt
import re
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.dates as mdates
import csv

class GeotailExtract(object):
    
    def __init__(self,filename,dti=None,dtf=None,autoplot=True,saveformat='png',mode='light',
                 x_visible=True,data_record=False):
        self.filename = filename #file name 
        self.dti = dti #start of datetime range
        self.dtf = dtf #end of datetime range
        self.saveformat = saveformat #Format of final image to be saved
        
        self.content,self.instr = self.gzread(self.filename) #Use gzread to open and extract data from file
        
        #Extract data again to make arrays
        self.dtime_arr,self.energy_arr,self.count_data = self.LEP_3Dextract(self.content)
        
        #Filter out data by datetime window
        if self.dti != None and self.dtf != None: #Currently only available for double-end filter i.e., there are dti and dtf.
            self.dtime_filter(self.dti, self.dtf)
        
        #Plot LEP Data
        if autoplot == True:
            self.LEP_3DPlotter(self.dtime_arr,self.energy_arr,self.count_data,self.saveformat,mode,x_visible,data_record)
    
    def gzread(self,filename):
        #gzip file unpacking function. Applicable only to LEP file
        
        if 'lep' in filename:
            instr = 'LEP'
        
        file = gzip.open(filename,'rb')
        content = file.read()
        
        return content, instr
    
    def LEP_3Dextract(self,content):
        #Method to read and obtain data from opened LEP data file (Editor A)

        #Split file content into lines
        content_lines = content.splitlines()
        
        #Total number of lines
        no_lines = len(content_lines)
        
        #List to append
        dtime_list = []
        
        #Dictionaries for energy and count
        dict1 = dict() #Reset with every datetime
        dict2 = dict() #Final with format key: [sum for all datetimes]
        
        for i in np.arange(no_lines):
            #Length of line can be used to separate which type of data it contains
            
            if len(content_lines[i]) == 17: #Datetime
                #Cannot use indexing for bytes, only slicing
                dtime_str = str(content_lines[i]) #Convert to string
                dtime_str = dtime_str[2:19] #Strip out the \b and final '
                dtime = dt.datetime.strptime(dtime_str,'%Y%m%d %H:%M:%S') #Make dt.datetime object
                dtime_list.append(dtime) #Append to list
            elif len(content_lines[i]) == 64: #Data
                #Get energy and count of this line
                data_line = str(content_lines[i])
                data_line = re.sub(' +',' ',data_line) #Make all columns separated by one space
                data_split = data_line.split(' ') #Split data, first element is 'b' for byte line break
                energy = data_split[1] #Keep string for dictionary keys
                count = float(data_split[5])
                
                #Append data
                if energy not in dict1.keys():
                    if i < 3586: #The very 1st datetime, format is different
                        dict1.update({energy:[count]}) #Make list with count as 1st element
                    else: #For next datetimes, should have empty list ready to append. Assume energy keys do not change.
                        dict1[energy].append(count)
                else:
                    dict1[energy].append(count) #Append existing list
                
            elif len(content_lines[i]) == 0: #If there is no data.
                #Calculate count for each energy and update
                for key in dict1.keys():
                    total = np.sum(dict1[key]) #Can sum list
                    #ave = total/16
                    if key not in dict2.keys():
                        dict2.update({key:[total]})
                    else:
                        dict2[key].append(total)
                    
                    #Reset dictionary arrays, but not keys
                    dict1[key] = []

        #Make sure datetime list is array
        dtime_arr = np.array(dtime_list)
        
        #Obtain energy array
        energy_list = dict2.keys()
        energy_arr = np.array([float(e) for e in energy_list])
        
        #Obtain (energy, datetime) count array
        count_data = np.nan*np.zeros((len(energy_arr),len(dtime_arr)))
        for key in dict2.keys():
            count_arr = np.array(dict2[key])
            if len(count_arr) != len(dtime_arr):
                raise Exception('Count missing for some datetime')
            else:
                row = np.where(energy_arr == float(key))[0][0]
                count_data[row,:] = count_arr #Set whole row to count array
        
        return dtime_arr,energy_arr,count_data #dtime, keV, count
    
    def dtime_filter(self,dti,dtf):
        #Filter data from file by datetime range
        dt_copy = np.copy(self.dtime_arr) #Copy array to not mess up assignment
        energy_copy = np.copy(self.energy_arr)
        count_copy = np.copy(self.count_data)
        
        if dti < dt_copy[0]: #If datetime range starts before data availability, fill array with blank rows
            dt_prior = []
            s_iter = np.floor(abs(dt_copy[0]-dti).seconds/12) #Find number of blank rows to fill
            for n in np.arange(s_iter):
                dt_prior.append(dti + dt.timedelta(seconds=n*12)) #Create a list of blank datetimes
            
            dt_prior = np.array(dt_prior) #Make the list a NumPy array
            count_prior = np.nan*np.zeros((len(energy_copy),len(dt_prior))) #Fill count data array with NaNs for blank datetimes.
            dt_copy = np.concatenate((dt_prior,dt_copy)) #Connect blank datetime array with data datetime array.
            count_copy = np.concatenate((count_prior,count_copy),axis=1) #Connect blank data array with data array
        
        
        if dtf > dt_copy[-1]: #If datetime range finishes after data availability, fill array with blank rows (same method as above).
            dt_after = []
            s_iter = np.floor(abs(dtf-dt_copy[-1]).seconds/12)
            for n in np.arange(s_iter): #count down from dtf to keep limits fixed
                dt.after.append(dtf - dt.timedelta(seconds=n*12))
            
            dt_after = np.array(dt_after)
            count_after = np.nan*np.zeros((len(energy_copy),len(dt_after)))
            
            dt_copy = np.concatenate((dt_copy,dt_after))
            count_copy = np.concatenate((count_copy,count_after),axis=1)
            
        #Filter out datetimes outside of defined range
        dt_copy = dt_copy[dt_copy > dti]
        dt_copy = dt_copy[dt_copy < dtf]
        
        #Find row index in data array corresponding to the defined datetime range.
        low_ind = np.where(dt_copy > dti)[0][0] #Lower indexing Inclusive
        high_ind = np.where(dt_copy > dti)[0][-1] + 1 #Upper indexing exclusive
        
        #Filter energy and count data arrays
        energy_copy = energy_copy[low_ind:high_ind]
        count_copy = count_copy[:,low_ind:high_ind]
        
        #Update filtered arrays to self instances.
        self.dtime_arr = dt_copy
        self.energy_arr = energy_copy
        self.count_data = count_copy
        
    def comb_array(self,dtime_arr,energy_arr,count_data):
        #Create combined array in format of Cluster plotter (see: PanelPlotter)
        if len(dtime_arr) >= len(energy_arr):
            nan_extra = np.nan * \
                np.zeros((len(dtime_arr)-len(energy_arr)))
            energy_arr = np.concatenate((energy_arr, nan_extra))
        else:
            nan_extra = np.nan * \
                np.zeros((len(energy_arr)-len(dtime_arr)))
            nan_extra2D = np.nan * \
                np.zeros(
                    (len(energy_arr)-len(dtime_arr), len(count_data[0, :])))

            dtime_arr = np.concatenate((dtime_arr, nan_extra))
            count_data = np.concatenate((count_data, nan_extra2D), axis=0)
            
        comb1 = np.stack(
            (dtime_arr, energy_arr)).transpose()
        #Stitch data along column
        comb2 = np.concatenate(
            (comb1, count_data.transpose()), axis=1) #This is different from Cluster
        
        return comb2
    
    
    def LEP_3DPlotter(self,dtime_arr,energy_arr,count_data,saveformat,mode,x_visible,data_record=False):
        #Plotting function for Geotail's ion distribution flux
        
        #Set relevant fontsizes
        subtitle_fontsize = 32
        axlabel_fontsize = 24
        axtick_fontsize = 20
        cblabel_fontsize = 24
        cbtick_fontsize = 20

        #Set color map
        cmap_name = 'nipy_spectral'
        
        #Define figure properties for light and dark mode
        if mode == 'light':
            ctick = 'k' #Black ticks/labels
            tchoice = False  #Plot with white background
        elif mode == 'dark':
            ctick = 'w' #White ticks/labels
            tchoice = True #Plot with transparent background
        else: #Raise exception for incorrect input.
            raise Exception('Incorrect Colour Mode')
        
        #Inverse extremities of nipy_spectral for light mode only.
        if mode == 'light':
            cm = plt.get_cmap(cmap_name,10).copy() #Get nipy_spectral with quantisation level of 10. Copy instance to not impact global settings
            clist = [] #Create list of 10 color tuples
            for a in range(10):
                clist.append(cm(a)) #Append each color tuple (call using curve bracket on cm, not square)
            clist2 = clist.copy() #Copy list for new color map.
            clist2[0] = clist[-1] #Swap first and last (black and grey) color coordinates.
            clist2[-1] = clist[0]
            n_bins = 256 #Set quantisation level of new colour map, same as before.
            cmap_name2 = 'nipy_swapped' #Set new colormap name
            cm = colors.LinearSegmentedColormap.from_list(cmap_name2,clist2,n_bins) #Create new colormap and set to this name
        else:
            cm = plt.get_cmap(cmap_name) #In dark mode, keep using nipy_spectral.
        
        #Create grid mesh and plot
        h,v = np.meshgrid(dtime_arr,energy_arr)
        fig,ax = plt.subplots(1,1,figsize=(20,5))
        im = ax.pcolormesh(h,v,count_data,cmap=cm)#
            
        #Obtain x-axis labels and names (10 min per major tick, 1 min per minor tick, similar to PanelPlotter.py).
        dtime_ticks = [] 
        dtime_ticklabels = []
        dtime_range = (self.dtf - self.dti).seconds/600 #Use dtf and dti for defined datetime range.
        m=0
        while m <= dtime_range:
            tick_dt = dtime_arr[0] + m*dt.timedelta(seconds=600)
            dtime_ticks.append(tick_dt)
            tick_time_split = str(tick_dt).split(' ')[1].split(':')
            tick_dt_str = ':'.join((tick_time_split[0],tick_time_split[1]))
            dtime_ticklabels.append(tick_dt_str)
            m = m+1
        
        #x-axis label settings (Similar to PanelPlotter.py)
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(AutoMinorLocator(10))
        ax.tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
        ax.tick_params(axis='x',which='major',direction='inout',width=1,length=15)
        ax.axes.get_xaxis().set_ticks(dtime_ticks,labels=[])
        ax.tick_params(axis='x',which='minor',direction='inout',width=1,length=5)
        dt_str = dtime_arr[0].strftime('%Y-%m-%d') + ' (UT)'
        
        #Set x-axis labels visible only when dictated
        if x_visible == True:
            ax.axes.get_xaxis().set_ticks(dtime_ticks,labels=dtime_ticklabels)
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            ax.set_xlabel(dt_str,fontsize=axlabel_fontsize,fontfamily='sans-serif') 
        
        #y-axis label settings
        ax.set_yscale('log')
        ax.tick_params(axis='y',which='major',direction='inout',width=1,length=15) #if log, use default minor ticks location
        ax.tick_params(axis='y',which='minor',direction='inout',width=1,length=5)
        #ax.set_ylim(vert_min, vert_max)
        ax.set_ylabel('Energy (keV)',fontsize=axlabel_fontsize,fontfamily='sans-serif')
        
        #Color bar settings
        axins = inset_axes(ax,width="1%",height="100%",loc='right',borderpad=-2)
        cb = fig.colorbar(im,cax=axins,orientation='vertical',pad=0.01)
        cb.set_label(label='Count (unit)',size=cblabel_fontsize,fontfamily='sans-serif')
        cb.ax.tick_params(labelsize=cbtick_fontsize)
        ax.tick_params(axis='both', which='major', labelsize=axtick_fontsize)
        
        #In-plot label
        ax.text(0.04,0.85,'(F)',fontsize=subtitle_fontsize,fontfamily='sans-serif',
                   horizontalalignment='center',verticalalignment='center',
                   transform=ax.transAxes,bbox=dict(facecolor='white',edgecolor='None',alpha=0.8))
        ax.text(0.93, 0.85, 'Geotail', fontsize=subtitle_fontsize,fontfamily='sans-serif',
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
        
        #Set colors
        ax.spines['bottom'].set_color(ctick) #Set axes color
        ax.spines['top'].set_color(ctick)
        ax.spines['left'].set_color(ctick)
        ax.spines['right'].set_color(ctick)
        ax.tick_params(axis='both',which='both',colors=ctick) #Set tick color
        ax.xaxis.label.set_color(ctick) #Set x-axis ticklabel color
        ax.yaxis.label.set_color(ctick) #Set y-axis ticklabel color
        cb.ax.axes.tick_params(axis='y',which='both',colors=ctick) #Set colorbar tick color
        cb.ax.axes.yaxis.label.set_color(ctick) #Set colorbar label color
        
        #Save figure
        savename = 'Geotail_LEP_{}.{}'.format(dt_str.split(' ')[0],saveformat)
        plt.subplots_adjust(left=0.08,bottom=0.16,right=0.90,top=0.97,
                            wspace=0.0,hspace=0.08)
        plt.savefig(savename,dpi=300,format=saveformat,transparent=tchoice)
        plt.close()
        
        #Write to CSV
        if data_record == True:
            comb_arr = self.comb_array(dtime_arr,energy_arr,count_data)
            start_dtstr = dtime_arr[0].strftime('%Y-%m-%dT%H-%M-%SZ')
            stop_dtstr = dtime_arr[-1].strftime('%Y-%m-%dT%H-%M-%SZ')
            csv_filename = '_'.join(['IonFlux','GEOTAIL',start_dtstr,stop_dtstr]) + '.csv'
            #see CSA_InstrData for code
            with open(csv_filename,'w', newline = '') as csvfile:
                csvwriter = csv.writer(csvfile)
            
                for i in np.arange(len(comb_arr[:,0])):
                    csvwriter.writerow(comb_arr[i,:])
            csvfile.close()
            print('IonFlux, GEOTAIL data recorded in csv')
    
"""
Example call script
"""
"""
filename = '20020318_lep_psd_2628.txt.gz'
dti = dt.datetime(2002,3,18,14,15,0)
dtf = dt.datetime(2002,3,18,15,15,0)
GO = GeotailExtract(filename,dti,dtf,autoplot=True,mode='light',saveformat='png',
                    x_visible=True,data_record=True)
"""

