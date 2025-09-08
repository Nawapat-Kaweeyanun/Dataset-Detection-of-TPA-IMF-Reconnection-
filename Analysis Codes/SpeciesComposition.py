"""
Edition Date: 2025-August-29
@author: Nawapat Kaweeyanun
"""

"""
Objective: Script for O+/H+ ratio using Cluster CODIF data.
"""

import numpy as np
import csv
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import AutoMinorLocator
from string import ascii_uppercase as ALC
from PanelPlotter import CPanel

class SpeciesComp(object):

    def __init__(self,start_dt,stop_dt,sc,saveformat,pa_spec=None,new_dl=True,
                 autoplot=True,label_start='A',xax_view=0,mode='light',
                 data_record=False):
        
        #Method for computing and plotting the ratio between O+/H+ ions.
        #Note: parameter name has to be 'ObyH'.
        
        self.start_dt = start_dt #Period start datetime
        self.stop_dt = stop_dt #Period end datetime
        self.sc = sc #Spacecraft (only one - C4)
        self.saveformat = saveformat #Format of figure file to save
        self.pa_spec = pa_spec #Pitch angle, if specified
        
        #Forbid pa not 0-90-180 (as of current iteration). 
        if self.pa_spec not in [None,0,90,180]:
            raise Exception('Pitch angle not supported')
        
        #Check that figure format is png/svg
        if self.saveformat not in ['png','svg']:
            raise Exception('Figure format not supported')
        
        #Set up relevant speicies parameters
        self.paralist = ['H+_fluxU','He+_fluxU','O+_fluxU']
        
        #Select pitch angle if specified.
        if self.pa_spec != None:
            self.paralist = [para + str(pa_spec) for para in self.paralist]
        
        #Set up parameter dictionary (in format of PanelPlotter.py)
        self.paradict = {}
        self.paradict.update({sc:self.paralist})
        self.instrlist = ['CIS']
        
        #Call PanelPlotter up to Panel Arranger (no plotting)
        self.PPBundle = CPanel(self.start_dt,self.stop_dt,[self.sc],self.paralist,self.instrlist,self.saveformat,
                          new_dl,autoplot=False)
        
        #Obtain data dictionary from the CPanel object.
        self.DataDict = self.PPBundle.PanelDict
        print('Individual Parameters Data Obtained')

        #Call species ratio function.
        self.ObyHDict = self.ObyH(self.sc, self.paralist, self.pa_spec, self.DataDict)
        print('Fractional Arrays Obtained')
        
        #Find ALC index for first subplot
        self.sp_ind0 = ALC.find(label_start.upper())
        
        #Plot O+/H+ ratio data.
        if autoplot == True:
            self.CompPlotter(self.ObyHDict,self.sc,self.pa_spec,xax_view,self.saveformat,mode,data_record)
            print('Figure Plotted')
        
    def ObyH(self,sc,paralist,pa_spec,DataDict):
        #Method for computing the ratio between O+/H+ ions.
        
        #Obtain H+ and O+ data and isolate flux arrays.
        Hcomb_arr = DataDict['H+_fluxU'][sc]
        Ocomb_arr = DataDict['O+_fluxU'][sc]
        Hdt_arr,HE_arr,HF_arr = self.Extractor('H+_fluxU',Hcomb_arr)
        Odt_arr,OE_arr,OF_arr = self.Extractor('O+_fluxU',Ocomb_arr)
        
        #Note: H+ has 8-second resolution while O+ has 4-second resolution, but all H+ datetimes are in O+ datetimes (a subset).
        #Plot by H+ resolution to eliminate gaps in plots
        #Assume same energy tags for H+ and O+.  
        HF_arr2 = np.nan*np.zeros((len(HE_arr),len(Hdt_arr)))
        OF_arr2 = np.nan*np.zeros((len(HE_arr),len(Hdt_arr)))                
        
        for d in np.arange(len(Hdt_arr)):
            HF_dt = HF_arr[:,Hdt_arr==Hdt_arr[d]].ravel()
            HF_arr2[:,d] = HF_dt
            
            #H+ datetime is a subset of O, but millisecond differences exist.
            #Create closest match (works because O+ res is 4s and differences are less than 1s)
            diff_arr = abs(Odt_arr-Hdt_arr[d])
            Odt_ind = np.where(diff_arr == min(diff_arr))[0][0]
            OF_dt = OF_arr[:,Odt_ind].ravel()
            OF_arr2[:,d] = OF_dt
            
        
        #Divide O+ by H+ fluxes. Get inf (NaN when plotted) where H+ flux is zero.
        ObyH = OF_arr2/HF_arr2
        
        #Concatenate energy or datetime/count arrays to produce 2D combined arrays.
        if len(Hdt_arr) >= len(HE_arr):
            nan_extra = np.nan*np.zeros((len(Hdt_arr)-len(HE_arr)))
            HE_arr = np.concatenate((HE_arr,nan_extra))
        else:
            nan_extra = np.nan*np.zeros((len(HE_arr)-len(Hdt_arr)))
            nan_extra2D = np.nan*np.zeros((len(HE_arr)-len(Hdt_arr),len(Hdt_arr)))
            Hdt_arr = np.concatenate((Hdt_arr,nan_extra))
            ObyH = np.concatenate((ObyH,nan_extra2D),axis=0)
        
        comb1 = np.stack((Hdt_arr,HE_arr)).transpose() 
        comb2 = np.concatenate((comb1,ObyH.transpose()),axis=1) 
        PanelDict = {}
        if pa_spec != None:
            PanelDict.update({'ObyH{}'.format(pa_spec):comb2})
        else:
            PanelDict.update({'ObyH':comb2})
        return PanelDict
    
    
    def CompPlotter(self,PanelDict,sc,pa_spec,xax_view,saveformat,mode,data_record=False):
        #Method for plotting O+/H+ ratio data.
        
        #Obtain number of subplots and parameter list.
        no_subplots = len(PanelDict.keys())
        paralist = list(PanelDict.keys())
        
        
        #Set up figure and fontsizes.
        fig,ax = plt.subplots(no_subplots,1,figsize=(20,5*no_subplots),squeeze=False)
        subtitle_fontsize = 32
        axlabel_fontsize = 24
        axtick_fontsize = 20
        cblabel_fontsize = 18
        cbtick_fontsize = 18
        
        #Set up color mode.
        if mode == 'light':
            ctick = 'k'
            tchoice = False
        elif mode == 'dark':
            ctick = 'w'
            tchoice = True
        else:
            raise Exception('Incorrect Colour Mode')
        
        #Set up color map by light or dark mode (see PanelPlotter.py).
        cmap_name = 'nipy_spectral'
        if mode == 'light':
            cm = plt.get_cmap(cmap_name,10).copy() 
            clist = [] 
            for a in range(10):
                clist.append(cm(a))
            clist2 = clist.copy() 
            clist2[0] = clist[-1] 
            clist2[-1] = clist[0]
            n_bins = 256
            cmap_name2 = 'nipy_swapped' 
            cm = colors.LinearSegmentedColormap.from_list(cmap_name2,clist2,n_bins) 
        elif mode == 'dark':
            cm = plt.get_cmap(cmap_name)
        
        #Plot each panel
        for i in np.arange(no_subplots):
            para = paralist[i]
            comb_arr = PanelDict[para]
            dtime_arr,energy_arr,flux_arr = self.Extractor(para,comb_arr)
            h,v = np.meshgrid(dtime_arr,energy_arr)
            
            #Set up flux and energy limits and color bar labels
            flux_min = 0.0
            flux_max = 0.5 #anything over one will get plotted by irrelevant
            vert_min = 0.01
            vert_max = energy_arr.max()
            vert_label = 'Energy (keV)'
            if pa_spec == None:
                cblabel = 'O+/H+ Fraction'
            else:
                cblabel = 'O+/H+ Fraction (PA = ${}^\circ$)'.format(pa_spec)
            
            #Plot O+/H+ .
            im = ax[i,0].pcolormesh(h,v,flux_arr,vmin=flux_min,vmax=flux_max,cmap=cm)
                
            #Set up x-axis ticks (use initially defined datetime range).
            dtime_ticks = [] 
            dtime_ticklabels = []
            dtime_range = (self.stop_dt - self.start_dt).seconds/600 
            
            m=0
            while m <= dtime_range:
                tick_dt = self.start_dt + m*dt.timedelta(seconds=600)
                dtime_ticks.append(tick_dt)
                tick_time_split = str(tick_dt).split(' ')[1].split(':')
                tick_dt_str = ':'.join((tick_time_split[0],tick_time_split[1]))
                dtime_ticklabels.append(tick_dt_str)
                m = m+1
            
            #Set x-axis visibility by options (0 = bottommost subplot only, 1 = all invisible, 2 = all visible) and other properties.
            if xax_view == 0:
                if i < no_subplots-1:
                    ax[i,0].minorticks_on()
                    ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                    ax[i,0].tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=False)
                    ax[i,0].tick_params(axis='x',which='major',direction='inout',width=1,length=15)
                    ax[i,0].tick_params(axis='x',which='minor',direction='inout',width=1,length=5)
                    ax[i,0].set_xlim(self.start_dt,self.stop_dt) #Optional for plt.pcolormesh, but put here anyway
                else:
                    ax[i,0].minorticks_on()
                    ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                    ax[i,0].tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
                    ax[i,0].tick_params(axis='x',which='major',direction='inout',width=1,length=15)
                    ax[i,0].axes.get_xaxis().set_ticks(dtime_ticks,labels=dtime_ticklabels)
                    ax[i,0].tick_params(axis='x',which='minor',direction='inout',width=1,length=5)
                    ax[i,0].set_xlim(self.start_dt,self.stop_dt)
                    dt_str = self.start_dt.strftime('%Y-%m-%d') 
                    ax[i,0].set_xlabel('{} (UT)'.format(dt_str),fontsize=axlabel_fontsize)
            elif xax_view == 1:
                ax[i,0].minorticks_on()
                ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                ax[i,0].tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=False)
                ax[i,0].tick_params(axis='x',which='major',direction='inout',width=1,length=15)
                ax[i,0].tick_params(axis='x',which='minor',direction='inout',width=1,length=5)
                ax[i,0].axes.get_xaxis().set_ticks(dtime_ticks,labels=dtime_ticklabels)
                ax[i,0].set_xlim(self.start_dt,self.stop_dt)
            elif xax_view == 2:
                ax[i,0].minorticks_on()
                ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                ax[i,0].tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
                ax[i,0].tick_params(axis='x',which='major',direction='inout',width=1,length=15)
                ax[i,0].axes.get_xaxis().set_ticks(dtime_ticks,labels=dtime_ticklabels)
                ax[i,0].tick_params(axis='x',which='minor',direction='inout',width=1,length=5)
                ax[i,0].set_xlim(self.start_dt,self.stop_dt)
                dt_str = self.start_dt.strftime('%Y-%m-%d') 
                ax[i,0].set_xlabel('{} (UT)'.format(dt_str),fontsize=axlabel_fontsize)
                
            #Set y-axis properties.
            ax[i,0].set_yscale('log')
            ax[i,0].tick_params(axis='y',which='major',direction='inout',width=1,length=15) #if log, use default minor ticks location
            ax[i,0].tick_params(axis='y',which='minor',direction='inout',width=1,length=5)
            ax[i,0].set_ylim(vert_min, vert_max)
            ax[i,0].set_ylabel(vert_label,fontsize=axlabel_fontsize)
           
                
            #Set color bar properties and in-plot texts
            axins = inset_axes(ax[i,0],width="1%",height="100%",loc='right',borderpad=-2)
            cb = fig.colorbar(im,cax=axins,orientation='vertical',pad=0.01)
            cb.set_label(label=cblabel,size=cblabel_fontsize)
            cb.ax.tick_params(labelsize=cbtick_fontsize)
            ax[i,0].tick_params(axis='both', which='major', labelsize=axtick_fontsize)
            #Letter
            ax[i,0].text(0.04,0.85,'({})'.format(ALC[i+self.sp_ind0]),fontsize=subtitle_fontsize,
                       horizontalalignment='center',verticalalignment='center',
                       transform=ax[i,0].transAxes,bbox=dict(facecolor='white',edgecolor='None',alpha=0.5))
            #Spacecraft
            ax[i,0].text(0.96,0.85,sc,fontsize=subtitle_fontsize,
                       horizontalalignment='center',verticalalignment='center',
                       transform=ax[i,0].transAxes,bbox=dict(facecolor='white',edgecolor='None',alpha=0.5))
           
            #Set panel color properties.
            ax[i,0].spines['bottom'].set_color(ctick)
            ax[i,0].spines['top'].set_color(ctick)
            ax[i,0].spines['left'].set_color(ctick)
            ax[i,0].spines['right'].set_color(ctick)
            ax[i,0].tick_params(axis='both',which='both',colors=ctick)
            ax[i,0].xaxis.label.set_color(ctick)
            ax[i,0].yaxis.label.set_color(ctick) 
            cb.ax.axes.tick_params(axis='y',which='both',colors=ctick) 
            cb.ax.axes.yaxis.label.set_color(ctick)
            
            #Record data to CSV file.
            if data_record == True:
                self.csv_record(para, sc, comb_arr)
                print('{}, {} data recorded in csv'.format(para, sc))
        
        #Set up figure name.
        dt_str = str(self.start_dt).split(' ')[0]
        savename = 'SpeciesFraction_{}_{}'.format(dt_str,sc)
        for para in paralist:
            savename = savename + '_{}'.format(para)
        savename = savename + '.{}'.format(saveformat) 
        
        #Set subplot margins.
        plt.subplots_adjust(left=0.08,bottom=0.15,right=0.90,top=0.97,
                            wspace=0.0,hspace=0.08) #0.05 = 5% from left figure width/height
        
        #Save figure.
        plt.savefig(savename,dpi=96,format=saveformat,transparent=tchoice)
        plt.close()
  
        
    def Extractor(self,para,comb_arr):
        #Function to extract datetime/energy/data for each combined parameter array.
        
        #If there is no data, return blank lists.
        if np.size(comb_arr) == 0:
            print('No data for {}'.format(para))
            dtime_arr = []
            energy_arr = []
            fluxfloat_arr = []
            return dtime_arr,energy_arr,fluxfloat_arr
        
        #Remove NaN from datetime and energy arrays.
        dtime_arr = comb_arr[:,0]
        l1 = len(dtime_arr)
        dtime_arr = np.array([d for d in dtime_arr if type(d) is dt.datetime]) #drop nan elements (can't use isnan for dt)
        l2 = len(dtime_arr)
        energy_arr = comb_arr[:,1]
        energy_arr = np.array([e for e in energy_arr if np.isnan(e) == False]) #drop nan elements
        
        #Isolate flux data from combined array.    
        if l1 > l2:
            flux_arr = comb_arr[0:l2,2:].transpose()
        else:
            flux_arr = comb_arr[:,2:].transpose()
            
        #Convert flux into float arrays.    
        fluxfloat_arr = np.zeros(np.shape(flux_arr))
        for m in np.arange(len(flux_arr[:,0])):
            for n in np.arange(len(flux_arr[0,:])):
                fluxfloat_arr[m,n] = flux_arr[m,n]
                    
        return dtime_arr,energy_arr,fluxfloat_arr
    
    def csv_record(self,para,sc,comb_data):
        #Record data in a CSV file.
        start_dtstr = self.start_dt.strftime('%Y-%m-%dT%H-%M-%SZ')
        stop_dtstr = self.stop_dt.strftime('%Y-%m-%dT%H-%M-%SZ')
        csv_filename = '_'.join([para,sc,start_dtstr,stop_dtstr]) + '.csv'
        with open(csv_filename,'w', newline = '') as csvfile:
            csvwriter = csv.writer(csvfile)
                        
            for i in np.arange(len(comb_data[:,0])):
                csvwriter.writerow(comb_data[i,:])
                        
        csvfile.close()
    

"""
Example Script
"""
"""
start_dt = dt.datetime(2002,3,18,14,15,0)
stop_dt = dt.datetime(2002,3,18,15,15,0)
sc = 'C4'
saveformat = 'png'
pa_spec = None
Comp = SpeciesComp(start_dt, stop_dt, sc, saveformat,pa_spec,
                   new_dl=True,autoplot=True,label_start='D',xax_view=1,mode='light',data_record=True)
#xax_view = 0,1,2 (bottom-only,invisible,visible)
"""