"""
Edition Date: 2025-October-29
@author: Nawapat Kaweeyanun
"""

"""
Objective: Two classes are present in this file.
    (1) A data class for extracting low-resolution IMAGE data from a provided (ECLAT) data file.
    (2) A class for plotting high-resolution IMAGE data from a provided data file.

Prerequisite: cdflib module
"""

import os
import cdflib
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import datetime as dt
from scipy.io import readsav
import tarfile

"""
Class for extracting ECLAT IMAGE data file.
"""

class data(object):
    
    def __init__(self, cdffilename, DataID):
        
        #open file with cdflib
        cdffile = cdflib.CDF(cdffilename)
        
        #Extract parameters
        Epochdates = cdffile.varget('time_tags__' + DataID)
        self.datetimes = np.array(cdflib.cdfepoch.to_datetime(Epochdates))
        
        self.hemis = cdffile.varget('Hemisphere__' + DataID) #0=North, 1=South
        if type(self.hemis) == np.int32: #remove corrupted data
            self.hemis = np.array([np.nan])
        
        self.qual = cdffile.varget('Quality__' + DataID) #1=good, 0=potential inaccuracy
        
        self.pixelnum = cdffile.varget('Pixels__' + DataID)
        
        self.image = np.array(cdffile.varget('Image__' + DataID)) #10x40x40 XY array
        
        self.XD1 = np.array(cdffile.varget('dimension_x__' + DataID)) #40-size vector, in km
        
        self.YD2 = np.array(cdffile.varget('dimension_y__' + DataID)) #40-size vector, in km


"""
Create IMAGE-HighRes class object with targeted cdf file as input
"""

class ImHighRes(object):
    
    def __init__(self,idlfilelist,instr,dtlist):
        #Take .tar IMAGE file as input and call the extraction method.
        self.instr = instr
        self.DataDict = self.DataExtract(idlfilelist,self.instr,dtlist)
        
    def DataExtract(self,idlfilelist,instr,dtlist):
        #Extract data from IMAGE .tar file to create a data dictionary
        DataDict = {}
         
        #The provided file list contains data from more than one instruments.
        #Obtain files that match the instrument chosen.
        instr_filelist = []
        for fname in idlfilelist:
            if fname[0:3] == instr.lower(): #If first three letters of file name matches the instrument
                instr_filelist.append(fname)
        
        #Abort code if there is no file match (usually mean wrong file chosen or a typo in instrument name)
        if len(instr_filelist) == 0:
            raise Exception('No files match chosen instrument. Check inputs.')
        
        #Create datetime list corresponding to the file list (i.e., all datetimes in the file list)
        instr_dtarr = []
        for fname in instr_filelist:
            dtime_str = fname[3:14] #every file will have dtime in these indices
            dtime = dt.datetime.strptime(dtime_str,'%Y%j%H%M')
            instr_dtarr.append(dtime)
        instr_dtarr = np.array(instr_dtarr) #convert from list to array
        
        for dtp in dtlist:
            #For each datetime in the plotting list, find closest file match.
            
            #Compare plotting datetime to datetime list
            diff_arr = abs(dtp-instr_dtarr)
            min_ind = np.where(diff_arr == min(diff_arr))[0][0]
            dt_closest_str = instr_dtarr[min_ind].strftime('%Y%m%d%H%M%S')
            DataDict.update({dt_closest_str:{}})
            
            #Once closest file match is identified, read the file to get radiance, magnetic latitude, and MLT data.
            file_closest = instr_filelist[min_ind]
            filedata = readsav(file_closest)
            im = filedata['imageinfo']['image'][0]
            mlat = filedata['imageinfo']['mlat'][0]
            mlt = filedata['imageinfo']['mlt'][0]
            
            #Update relevant data to the data dictionary.
            DataDict[dt_closest_str].update({'image':im})
            DataDict[dt_closest_str].update({'mlat':mlat})
            DataDict[dt_closest_str].update({'mlt':mlt})
        
        return DataDict
    
    def Plotter(self,DataDict,timeseries=False,ncol=4,mode='light',cbar_on=True,
                save=True,savefolder='',saveformat='png',data_record=False):  
        #Function to make IMAGE plots from prescribed data dictionary.
        #Can plot all datetimes in a single figure (timeseries=True) or produce each datetime as a separate figure (timeseries=False)

        #Set axes color and figure transparency based on light or dark modes.
        if mode.lower() == 'light':
            ctick = 'k'
            tchoice = False
        elif mode.lower() == 'dark':
            ctick = 'w'
            tchoice = True
        else:
            raise Exception ('Incorrect Colour Mode')
            
        #Set color map limits according to instruments
        if self.instr == 'S12':
            vmin = 100
            vmax = 5000
            title_name = 'SI-12'
        elif self.instr == 'S13':
            vmin = 100
            vmax = 5000
            title_name = 'SI-13'
        elif self.instr == 'WIC':
            vmin = 500
            vmax = 5000
            title_name = 'WIC'
        
        #If plotting all datetimes as a time series.
        if timeseries == True:
            #Use datetimes from the data dictionary, each element is a string
            dtp_list = list(DataDict.keys())  
            no_subplots = len(dtp_list)
            
            #Find number of rows in the figure from the preset number of columns
            nrow = int(np.ceil(no_subplots/ncol)) 

            #Create figure whose width is defined by column number and whose height is defined by row number.
            fig,ax = plt.subplots(nrow,ncol,subplot_kw={'projection': 'polar'},figsize=[4*ncol,4*nrow],squeeze=False) 
            
            #Create array for identifying hemisphere for later overlay function
            hno_arr = []
            
            #Generate each subplot
            for i in np.arange(no_subplots):
                #Calculate row and column indices of the subplot.
                row_ind = int(np.floor(i/ncol))
                col_ind = np.remainder(i,ncol)

                #Due to upload order in the DataExtract function, the datetime list should be in chronological order
                #Retrieve data and unravel arrays to flat shapes.
                imp = DataDict[dtp_list[i]]['image'].ravel()
                mlatp = DataDict[dtp_list[i]]['mlat'].ravel()
                mlt = DataDict[dtp_list[i]]['mlt'].ravel()
                
                #Define hemisphere by year (2000-2003: N, 2004-2005: S). This is possible due to IMAGE's steady precession.
                if dtp_list[i][0:4] in ['2000','2001','2002','2003']:
                    hemis = 'N'
                    hno_arr.append(0)
                elif dtp_list[i][0:4] in ['2004','2005']:
                    hemis = 'S'
                    hno_arr.append(1)
                
                #Plot IMAGE data. Log-scale work better than linear scale here, but both options are available.
                sub_im = ax[row_ind,col_ind].scatter(np.pi+(mlt*np.pi/12),mlatp,c=imp,s=2,
                                                     norm=colors.LogNorm(vmin=vmin,vmax=vmax),
                                                     cmap='seismic',zorder=1)
                #sub_im = ax[row_ind,col_ind].scatter(np.pi+(mlt*np.pi/12),mlatp,c=imp,s=1,
                #                                     vmin=vmin,vmax=vmax,
                #                                     cmap='seismic',zorder=1)
                
                #Set angular axis properties
                ax[row_ind,col_ind].set_theta_zero_location('N')
                ax[row_ind,col_ind].set_rlim(90,40,5)
                xticks_loc = ax[row_ind,col_ind].get_xticks().tolist()
                ax[row_ind,col_ind].xaxis.set_major_locator(mticker.FixedLocator(xticks_loc))
                ax[row_ind,col_ind].set_xticklabels(['12','15','18','21','0','3','6','9'],fontsize=12)
                ax[row_ind,col_ind].spines['polar'].set_color(ctick)
                
                #Set radial axis properties and subplot title
                yticks_loc = ax[row_ind,col_ind].get_yticks().tolist()
                if hemis == 'N': #Labels depend on hemisphere.
                    ax[row_ind,col_ind].yaxis.set_major_locator(mticker.FixedLocator(yticks_loc))
                    ax[row_ind,col_ind].set_yticklabels(['','50','60','70','',''],fontsize=10)
                    dtp_formatted = dtp_list[i][8:10] + ':' + dtp_list[i][10:12]
                    ax[row_ind,col_ind].set_title('IMAGE {} North\n{}'.format(title_name,dtp_formatted) ,color=ctick,fontsize=12,fontfamily='sans-serif')
                elif hemis == 'S':
                    ax[row_ind,col_ind].yaxis.set_major_locator(mticker.FixedLocator(yticks_loc))
                    ax[row_ind,col_ind].set_yticklabels(['','-50','-60','-70','',''],fontsize=10)
                    dtp_formatted = dtp_list[i][8:10] + ':' + dtp_list[i][10:12]
                    ax[row_ind,col_ind].set_title('IMAGE {} South\n{}'.format(title_name,dtp_formatted) ,color=ctick,fontsize=12,fontfamily='sans-serif')
                ax[row_ind,col_ind].yaxis.label.set_color(ctick)
                ax[row_ind,col_ind].tick_params(axis='both',which='both',colors=ctick)
                ax[row_ind,col_ind].grid(axis='both',which='both',color=ctick)
            
                #Record subplot data to a CSV file if required.
                if data_record == True:
                    #Loop over each parameter
                    for para in ['image','mlat','mlt']:
                        data = DataDict[dtp_list[i]][para]
                        csv_filename = 'IMAGE_HighRes_Data_{}_{}_{}.csv'.format(self.instr,dtp_list[i],para)
                        with open(csv_filename,'w', newline = '') as csvfile:
                            csvwriter = csv.writer(csvfile)
                            for m in np.arange(len(data[:,0])):
                                csvwriter.writerow(data[m,:])
                        csvfile.close()
            
            #Set color bar properties (if a color bar is included) and figure layout.
            if cbar_on == True:
                fig.subplots_adjust(left=0.02,right=0.92,bottom=0.05,top=0.83,wspace=0.25,hspace=0.0)
                cbar_ax = fig.add_axes([0.94,0.10,0.02,0.75]) #set axis location
                cb = plt.colorbar(sub_im,cax=cbar_ax) #use colour range from final plots. All subplots have same limits.
                cb.set_label(label='Radiance (R)',size=12,family='sans-serif')
            else:
                plt.tight_layout()
            
            #Save figure
            if save == True:
                savename = 'IMAGE_{}_{}_{}_{}.{}'.format(self.instr,dtp_list[0],dtp_list[-1],hemis,saveformat)
                savedir = savefolder + '/{}/{}/{}/'.format(dtp_list[0][0:4],dtp_list[0][4:6],dtp_list[0][6:8])

                os.makedirs(savedir,exist_ok=True)
                savepath = savedir + savename
                plt.savefig(savepath,dpi=100,format=saveformat,transparent=tchoice)
                plt.show()
                plt.close(fig)
                return fig,ax,hno_arr

            elif save == False:
                plt.close(fig)
                return fig,ax,hno_arr
            
        #If plotting each datetime as a separate figure.    
        elif timeseries == False:
            
            for dtime_str in DataDict.keys():
                dtime_plot = dt.datetime.strptime(dtime_str,'%Y%m%d%H%M%S')
                
                #Create simple 4x4 in figure with only one subplot
                fig,ax = plt.subplots(1,1,subplot_kw={'projection': 'polar'},figsize=[4,4],squeeze=False)
                
                #Create array for identifying hemisphere for later overlay function
                hno_arr = []
                
                #Retrieve data and unravel arrays to flat shapes.
                imp = DataDict[dtime_str]['image'].ravel()
                mlatp = DataDict[dtime_str]['mlat'].ravel()
                mlt = DataDict[dtime_str]['mlt'].ravel()
                
                #Define hemisphere by year (2000-2003: N, 2004-2005: S). This is possible due to IMAGE's steady precession.
                if dtime_str[0:4] in ['2000','2001','2002','2003']:
                    hemis = 'N'
                    hno_arr.append(0)
                elif dtime_str[0:4] in ['2004','2005']:
                    hemis = 'S'
                    hno_arr.append(1)
        
                #Plot IMAGE data. Log-scale work better than linear scale here, but both options are available.
                sub_im = ax[0,0].scatter(np.pi+(mlt*np.pi/12),mlatp,c=imp,s=2,norm=colors.LogNorm(vmin=vmin,vmax=vmax),
                                cmap='seismic',zorder=1)
                #sub_im = ax[0,0].scatter(np.pi+(mlt*np.pi/12),mlatp,c=imp,s=1,vmin=vmin,vmax=vmax,cmap='seismic',zorder=1)
                
                #Set angular axis properties
                ax[0,0].set_theta_zero_location('N')
                ax[0,0].set_rlim(90,40,5)
                xticks_loc = ax[0,0].get_xticks().tolist()
                ax[0,0].xaxis.set_major_locator(mticker.FixedLocator(xticks_loc))
                ax[0,0].set_xticklabels(['12','15','18','21','0','3','6','9'],fontsize=12)
                ax[0,0].spines['polar'].set_color(ctick) #set angular axis color
                
                
                #Set radial axis properties and subplot title.
                yticks_loc = ax[0,0].get_yticks().tolist()
                if dtime_plot.year in [2000,2001,2002,2003]:
                    hemis = 0
                    ax[0,0].yaxis.set_major_locator(mticker.FixedLocator(yticks_loc))
                    ax[0,0].set_yticklabels(['','50','60','70','',''],fontsize=10)
                    ax[0,0].set_title('IMAGE {} North\n{}'.format(title_name,dtime_plot.strftime('%H:%M')),color=ctick,fontsize=12,fontfamily='sans-serif')
                elif dtime_plot.year in [2004,2005]:
                    hemis = 1
                    ax[0,0].yaxis.set_major_locator(mticker.FixedLocator(yticks_loc))
                    ax[0,0].set_yticklabels(['','-50','-60','-70','',''],fontsize=10)
                    ax[0,0].set_title('IMAGE {} South\n{}'.format(title_name,dtime_plot.strftime('%H:%M')),color=ctick,fontsize=12,fontfamily='sans-serif')
                ax[0,0].yaxis.label.set_color(ctick) #set y-axis ticklabel color
                ax[0,0].tick_params(axis='both',which='both',colors=ctick) #set tick color both axes
                ax[0,0].grid(axis='both',which='both',color=ctick) #set grid color
            
            
                #Set color bar properties (if a color bar is included) and figure layout.
                if cbar_on == True:
                    fig.subplots_adjust(left=0.10,right=0.75,bottom=0.05,top=0.83,wspace=0.0,hspace=0.0)
                    cbar_ax = fig.add_axes([0.84,0.10,0.02,0.75]) #set axis location
                    cb = plt.colorbar(sub_im,cax=cbar_ax) #use colour range from final plots. All subplots have same limits.
                    cb.set_label(label='Radiance (R)',size=12,family='sans-serif')
                else:
                    plt.tight_layout()
                    
                #Record subplot data to a CSV file if required.
                if data_record == True:
                    #Loop over each parameter
                    for para in ['image','mlat','mlt']:
                        data = DataDict[dtime_str][para]
                        csv_filename = 'IMAGE_HighRes_Data_{}_{}_{}.csv'.format(self.instr,dtime_str,para)
                        with open(csv_filename,'w', newline = '') as csvfile:
                            csvwriter = csv.writer(csvfile)
                            for i in np.arange(len(data[:,0])):
                                csvwriter.writerow(data[i,:])
                        csvfile.close()
                
                #Save figure
                if save == True:
                    savename = 'IMAGE_{}_{}_{}.{}'.format(self.instr,dtime_str,hemis,saveformat)
                    savedir = savefolder + '/{}/{}/{}/'.format(dtime_str[0:4],dtime_str[4:6],dtime_str[6:8])

                    os.makedirs(savedir,exist_ok=True)
                    savepath = savedir + savename
                    plt.savefig(savepath,dpi=100,format=saveformat,transparent=tchoice)
                    plt.show()
                    plt.close(fig)

                elif save == False:
                    plt.close(fig)

            #Return only final fig and ax as dummy variables.
            return fig,ax,hno_arr
