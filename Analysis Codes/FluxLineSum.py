"""
Edition Date: 2025-October-21
@author: Nawapat Kaweeyanun
"""

"""
Objective: Sum up top energy channels of input panels and display them as line series in a single panel.

Prerequisite: PanelPlotter.py module in the same directory.
"""

import os
import csv
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from string import ascii_uppercase as ALC
import PanelPlotter as PP

"""
Class Script
"""

class FluxLineSum(object):
    
    def __init__(self,dti,dtf,paralist,sc,energy_min,new_dl=True,just_sum=False,autoplot=True,label_start='A',
                 xax_view=0,mode='light',save=True,saveformat='png',data_record=False):
        self.dti = dti #Starting datetime
        self.dtf = dtf #Final datetime
        self.paralist = paralist #Parameter list where each is summed
        self.sc = sc #One Cluster spacecraft only
        self.energy_min = energy_min #Energy threshold above which fluxes are summed.
        
        #Only applied for flux-related parameters (with spacecraft potential as exception for photoelectron purposes)
        #Remove all other parameters from the list.
        for para in self.paralist:
            if 'flux' not in para and para != 'SC_Pot':
                self.paralist.remove(para)
                print('Parameter {} removed from list. Not flux related.'.format(para))
        
        #Set upper limit of parameters to sum to six. Threshold is 7 here because of spacecraft potential.
        if len(self.paralist) > 7:
            raise Exception('More than six parameters to sum. Reduce number of parameters')
        
        #Define parameters to call PanelPlotter.py
        self.paradict = {}
        self.paradict.update({sc:self.paralist})
        self.instrlist = ['CIS','PEA','EFW']
        
        #Call PanelPlotter up to PanelArranger and create dictionary of data.
        self.PPBundle = PP.CPanel(self.dti,self.dtf,[self.sc],self.paralist,self.instrlist,saveformat,
                          new_dl,autoplot=False)
        self.OGPanelDict = self.PPBundle.PanelDict #This will serve as data for future PanelDict.
        
        #Call method to add a 'LineSum' panel at the end of DataDict, or just the 'LineSum' panel itself
        #Use 'just_sum' toggle
        self.PanelDict = self.SumCreator(self.OGPanelDict,just_sum,self.energy_min)
        
        #Find ALC index for first subplot
        self.sp_ind0 = ALC.find(label_start.upper())
        
        #Call the plotting function if selected.
        if autoplot == True:
            self.LineSumPlotter(self.PanelDict,xax_view,mode,save,saveformat,data_record)
        
    def SumCreator(self,OGPanelDict,just_sum,energy_min):
        #Method to generate panel containing summed up parameters
        #PanelDict has 'LineSum' key entry that split into further dictionary with parameter keys.
        #This is so that all of 'LineSum' is in one panel, which means less complicated looping later.
        #Moreover, this means that each 'LineSum' line can have its own datetime and avoid any resolution issues.
        #Minimum energy threshold is set in keV.
        
        #Set up panel dictionary to include the original parameters or not.
        if just_sum == False: #Include original parameters
            PanelDict = OGPanelDict
            PanelDict.update({'LineSum':{self.sc:{}}}) #Update dictionary with 'LineSum' key
        elif just_sum == True: #Exclude original parameters.
            PanelDict = {}
            PanelDict.update({'LineSum':{self.sc:{}}})
        
        #Retrieve data using the Extractor method
        for para in self.paralist:
            comb_arr = OGPanelDict[para][self.sc]
            dtime_arr,energy_arr,flux_arr = self.Extractor(para,comb_arr)
            
            #Energy in descending order, so earliest columns belongs to highest energy
            #Filter out flux so that only columns above the threshold remain
            Ebool = energy_arr >= energy_min
            flux_filtered = flux_arr[Ebool,:] #The Extractor makes datetime row, not column
            
            #Sum flux across rows
            flux_sum = np.sum(flux_filtered,axis=0) #has same rows as dtime_arr
            
            #Stack with datetime to form plotting array
            plot_arr = np.vstack((dtime_arr,flux_sum)).transpose()
            
            #Update to PanelDict['LineSum'] with 'LineSum_Para' as key
            PanelDict['LineSum'][self.sc].update({para:plot_arr})
        
        return PanelDict
    
    def Extractor(self,para,comb_arr):
        #Extract datetime, energy, and data for each data panel.    
        if np.size(comb_arr) == 0:
            print('No data for {}'.format(para))
            dtime_arr = []
            energy_arr = []
            fluxfloat_arr = []
            return dtime_arr,energy_arr,fluxfloat_arr
        
        #Remove NaN elements from datetime and energy arrays.
        dtime_arr = comb_arr[:,0]
        l1 = len(dtime_arr)
        dtime_arr = np.array([d for d in dtime_arr if type(d) is dt.datetime]) 
        l2 = len(dtime_arr)
        energy_arr = comb_arr[:,1]
        energy_arr = np.array([e for e in energy_arr if np.isnan(e) == False])
        if l1 > l2: #If there has been reduction in datetime array lengths, removed NaNs stacked at the end.
            flux_arr = comb_arr[0:l2,2:].transpose()
        else:
            flux_arr = comb_arr[:,2:].transpose()
            
        #Convert numbers in flux arrays to floats.    
        fluxfloat_arr = np.zeros(np.shape(flux_arr))
        for m in np.arange(len(flux_arr[:,0])):
            for n in np.arange(len(flux_arr[0,:])):
                fluxfloat_arr[m,n] = flux_arr[m,n]
                    
        return dtime_arr,energy_arr,fluxfloat_arr
    
    def LineSumPlotter(self,PanelDict,xax_view,mode,save,saveformat,data_record):
        #Plotting function. Same as Plot_Func of PanelPlotter.py except with extra LineSum functionality.
        
        #Obtain number of subplots and parameter lists.
        no_subplots = len(PanelDict.keys())
        paralist = list(PanelDict.keys()) #No longer the same list as the initial self.paralist
        
        #Subtract one panel if sc potential is used. (Cannot subtract from paralist itself because needs the data earlier.)
        if 'SC_Pot' in paralist:
            no_subplots = no_subplots - 1
        
        #Set up figure
        fig,ax = plt.subplots(no_subplots,1,figsize=(20,5*no_subplots),squeeze=False)
        subtitle_fontsize = 32
        axlabel_fontsize = 24
        axtick_fontsize = 20
        cblabel_fontsize = 18
        cbtick_fontsize = 18
        legend_fontsize = 14
        
        #Mode color settings
        if mode == 'light':
            ctick = 'k'
            tchoice = False
        elif mode == 'dark':
            ctick = 'w'
            tchoice = True
        else:
            raise Exception('Incorrect Colour Mode')
            
        #Mode color map settings (reverse if dark mode)
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
            
            #For original parameters, same principles as PanelPlotter.py.
            if 'LineSum' not in para:
                #Extract data using the Extractor Method
                comb_arr = PanelDict[para][self.sc]
                dtime_arr,energy_arr,flux_arr = self.Extractor(para,comb_arr)
                h,v = np.meshgrid(dtime_arr,energy_arr)
                
                #Set limits and labels.
                if para in ['ion_fluxU']:
                    flux_min = 10**4 
                    flux_max = 10**9 
                    vert_min = 0.01
                    vert_max = energy_arr.max()
                    vert_label = 'Energy (keV)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['H+_fluxU', 'He+_fluxU', 'O+_fluxU']:
                    flux_min = 10**3 
                    flux_max = 10**8
                    vert_min = 0.01
                    vert_max = energy_arr.max()
                    vert_label = 'Energy (keV)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['ion_fluxU0', 'ion_fluxU90', 'ion_fluxU180']:
                    flux_min = 10**6
                    flux_max = 10**9
                    vert_min = 5
                    vert_max = energy_arr.max()
                    vert_label = 'Energy (keV)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF (PA = {})'.format(para.split('U')[-1]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['H+_fluxU0', 'H+_fluxU90', 'H+_fluxU180',
                              'He+_fluxU0', 'He+_fluxU90', 'He+_fluxU180',
                              'O+_fluxU0', 'O+_fluxU90', 'O+_fluxU180']:
                    flux_min = 10**3 
                    flux_max = 10**6 
                    vert_min = 0.01
                    vert_max = energy_arr.max()
                    vert_label = 'Energy (keV)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF (PA = {})'.format(para.split('U')[-1]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['ion_fluxPA', 'H+_fluxPA', 'O+_fluxPA', 'He+_fluxPA']:
                    flux_min = 2*10**3 
                    flux_max = 2*10**7 
                    vert_min = 0
                    vert_max = 180
                    vert_label = 'Pitch Angle ($^\circ$)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['ion_fluxU_Sun', 'ion_fluxU_Dawn', 'ion_fluxU_Antisun', 'ion_fluxU_Dusk']:
                    flux_min = 10**4
                    flux_max = 10**5
                    vert_min = 0.01
                    vert_max = energy_arr.max()
                    vert_label = 'Energy (keV)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF ({}ward)'.format(para.split('_')[2]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk',
                              'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk',
                              'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk']:
                    flux_min = 10**4
                    flux_max = 10**6
                    vert_min = 0.05
                    vert_max = energy_arr.max()
                    vert_label = 'Energy (keV)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF ({}ward)'.format(para.split('_')[2]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['e_fluxU']:
                    flux_min = 5*10**6 #5*10**6
                    flux_max = 5*10**9
                    vert_min = 5/1000
                    vert_max = energy_arr.max()
                    vert_label = 'Energy (keV)'
                    cblabel = 'Electron DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['e_fluxU0', 'e_fluxU90', 'e_fluxU180']:
                    flux_min = 5*10**6
                    flux_max = 3*10**9
                    vert_min = 5/1000
                    vert_max = energy_arr.max()
                    vert_label = 'Energy (keV)'
                    cblabel = 'Electron DEF (PA = {})'.format(para.split('U')[-1]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para == 'e_fluxPA':
                    flux_min = 10**6 
                    flux_max = 10**8 
                    vert_min = 0
                    vert_max = 180
                    vert_label = 'Pitch Angle ($^\circ$)'
                    cblabel = 'Electron DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                
                #Plot panel data
                im = ax[i,0].pcolormesh(h,v,flux_arr,norm=colors.LogNorm(vmin=flux_min,vmax=flux_max),cmap=cm,edgecolors='none')
                
                #Plot photoelectron thresholds and spacecraft potentials.
                if 'SC_Pot' in paralist and 'e_fluxU' in para:
                    photoe_dtime = PanelDict['SC_Pot'][self.sc][:, 0]
                    photoe_energy = -1 * \
                        PanelDict['SC_Pot'][self.sc][:, 1]  
                    photoe_energy = photoe_energy/1000 
    
                    ax[i,0].plot(photoe_dtime, photoe_energy,
                                 c='k',ls='solid', linewidth=2)

                    ax[i,0].plot(photoe_dtime, photoe_energy *
                                 1.30, c='k',ls=(5,(10,5)), linewidth=2)

                
                #Set up x-axis ticks and labels
                dtime_ticks = []  
                dtime_ticklabels = []
                if (self.dtf-self.dti).seconds > 1200:
                    major_tick_interval = 600
                else:
                    major_tick_interval = 60
                dtime_range = (self.dtf - self.dti).seconds/major_tick_interval
                m = 0
                while m <= dtime_range:
                    tick_dt = self.dti + m*dt.timedelta(seconds=major_tick_interval)
                    dtime_ticks.append(tick_dt)
                    tick_time_split = str(tick_dt).split(' ')[1].split(':')
                    tick_dt_str = ':'.join(
                        (tick_time_split[0], tick_time_split[1]))
                    dtime_ticklabels.append(tick_dt_str)
                    m = m+1

                #Set x-axis properties.
                if xax_view == 0: #Label visible only in bottom-most plots.
                    if i < no_subplots-1:
                        ax[i,0].minorticks_on()
                        ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                        ax[i,0].tick_params(
                            axis='x', which='both', bottom=True, top=False, labelbottom=False)
                        ax[i,0].tick_params(
                            axis='x', which='major', direction='inout', width=1, length=15)
                        ax[i,0].axes.get_xaxis().set_ticks(
                            dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize) 
                        ax[i,0].tick_params(
                            axis='x', which='minor', direction='inout', width=1, length=5)
                        ax[i,0].set_xlim(self.dti, self.dtf)
                    else:
                        ax[i,0].minorticks_on()
                        ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                        ax[i,0].tick_params(
                            axis='x', which='both', bottom=True, top=False, labelbottom=True)
                        ax[i,0].tick_params(
                            axis='x', which='major', direction='inout', width=1, length=15)
                        ax[i,0].axes.get_xaxis().set_ticks(
                            dtime_ticks, labels=dtime_ticklabels)
                        ax[i,0].tick_params(
                            axis='x', which='minor', direction='inout', width=1, length=5)
                        ax[i,0].set_xlim(self.dti, self.dtf)
                        dt_str = self.dti.strftime(
                            '%Y-%m-%d')  # set date string
                        ax[i,0].set_xlabel('{} (UT)'.format(
                            dt_str), fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif xax_view == 1: #All x-axis labels invisible.
                    ax[i,0].minorticks_on()
                    ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                    ax[i,0].tick_params(
                        axis='x', which='both', bottom=True, top=False, labelbottom=False)
                    ax[i,0].tick_params(
                        axis='x', which='major', direction='inout', width=1, length=15)
                    ax[i,0].axes.get_xaxis().set_ticks(
                        dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize) #must set to ensure ticks are in correct position
                    ax[i,0].tick_params(
                        axis='x', which='minor', direction='inout', width=1, length=5)
                    ax[i,0].set_xlim(self.dti, self.dtf)
                elif xax_view == 2: #All x-axis labels visible.
                    ax[i,0].minorticks_on()
                    ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                    ax[i,0].tick_params(
                        axis='x', which='both', bottom=True, top=False, labelbottom=True)
                    ax[i,0].tick_params(
                        axis='x', which='major', direction='inout', width=1, length=15)
                    ax[i,0].axes.get_xaxis().set_ticks(
                        dtime_ticks, labels=dtime_ticklabels)
                    ax[i,0].tick_params(
                        axis='x', which='minor', direction='inout', width=1, length=5)
                    ax[i,0].set_xlim(self.dti, self.dtf)
                    dt_str = self.dti.strftime(
                        '%Y-%m-%d')  # set date string
                    ax[i,0].set_xlabel('{} (UT)'.format(
                        dt_str), fontsize=axlabel_fontsize, fontfamily='sans-serif')
                    
                #Set y-axis properties.
                if 'U' in para:
                    ax[i,0].tick_params(
                        axis='y', which='major', direction='inout', width=1, length=15)
                    ax[i,0].tick_params(
                        axis='y', which='minor', direction='inout', width=1, length=5)
                    ax[i,0].set_ylabel(
                        vert_label, fontsize=axlabel_fontsize, fontfamily='sans-serif')
                    ax[i,0].set_ylim(vert_min, vert_max)
                    ax[i,0].set_yscale('log')
                elif 'PA' in para:
                    ax[i,0].yaxis.set_minor_locator(AutoMinorLocator(5))
                    ax[i,0].tick_params(
                        axis='y', which='major', direction='inout', width=1, length=15)
                    ax[i,0].tick_params(
                        axis='y', which='minor', direction='inout', width=1, length=5)
                    ax[i,0].set_ylim(vert_min, vert_max)
                    ax[i,0].set_ylabel(
                        vert_label, fontsize=axlabel_fontsize, fontfamily='sans-serif')
                
                #Set color bar location and properties.
                axins = inset_axes(
                    ax[i,0], width="1%", height="100%", loc='right', borderpad=-2)
                cb = fig.colorbar(
                    im, cax=axins, orientation='vertical', pad=0.01)
                cb.set_label(label=cblabel, size=cblabel_fontsize, fontfamily='sans-serif')
                cb.ax.tick_params(labelsize=cbtick_fontsize)
                ax[i,0].tick_params(
                    axis='both', which='major', labelsize=axtick_fontsize)
                
                #Set in-plot text labels
                #Panel letters
                ax[i,0].text(0.04, 0.85, '({})'.format(ALC[i+self.sp_ind0]), fontsize=subtitle_fontsize, fontfamily='sans-serif',
                             horizontalalignment='center', verticalalignment='center',
                             transform=ax[i,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                #Spacecraft indicators
                ax[i,0].text(0.96, 0.85, self.sc, fontsize=subtitle_fontsize, fontfamily='sans-serif',
                             horizontalalignment='center', verticalalignment='center',
                             transform=ax[i,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                    
                #Set axes colors depending on light or dark modes.
                ax[i,0].spines['bottom'].set_color(ctick)
                ax[i,0].spines['top'].set_color(ctick)
                ax[i,0].spines['left'].set_color(ctick)
                ax[i,0].spines['right'].set_color(ctick)
                ax[i,0].tick_params(axis='both',which='both',colors=ctick)
                ax[i,0].xaxis.label.set_color(ctick) 
                ax[i,0].yaxis.label.set_color(ctick) 
                cb.ax.axes.tick_params(axis='y',which='both',colors=ctick) 
                cb.ax.axes.yaxis.label.set_color(ctick) 
                
                #Record data to CSV files if required.
                if data_record == True:
                    self.csv_record(para, self.sc, comb_arr)
                    print('{}, {} data recorded in csv'.format(para, self.sc))
                    
            #Plot sum of fluxes above certain energy threshold for each parameter.
            elif para == 'LineSum':
                
                #Set line colors (up to six)
                color_list = ['black','red','blue','green','magenta','gold']
                

                c = 0 #Initial color index for the color list. Add 1 per each loop over parameters below.
                
                #Loop over parameters (keys in each LineSum minidictionary)
                for key in PanelDict['LineSum'][self.sc].keys():
                    
                    #Obtain datetime and flux data for each parameter.
                    dt_arr = PanelDict['LineSum'][self.sc][key][:,0]
                    data_arr = PanelDict['LineSum'][self.sc][key][:,1]
                    
                    #List of parameters that use 3D flux distributions rather than pitch-angle flux distributions (PAD).
                    Para3D = ['H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk',
                                  'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk',
                                  'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk']
                    
                    #Set name list for legend for each type of parameters
                    if 'fluxU' in key and key not in Para3D: #PAD flux sorted by energy
                        if key.split('U')[-1] in ['0','90','180']: #If by-energy flux is at specific pitch angles.
                            paraname = key.split('_')[0] + ' PA = {}'.format(key.split('U')[-1]) 
                        else:
                            paraname = key.split('_')[0]
                    elif 'fluxU' in key and key in Para3D: #3D flux sorted by energy
                        paraname = key.split('_')[0] + ' {}'.format(key.split('_')[-1])
                    elif 'fluxPA' in key: #Flux sorted by pitch angles.
                        paraname = key.split('_')[0] + ' (PA)'
                    
                    #Plot summed flux line.
                    ax[i,0].plot(dt_arr, data_arr, label=paraname,
                                 color=color_list[c], ls='-', lw=2, alpha=1, zorder=2)
                    c = c+1 #Add one to the color index for next time.
                    
                    #Add a line at zero if needed.
                    if data_arr.min() < 0 and data_arr.max() > 0:
                        ax[i,0].axhline(y=0,color='k',lw=1,alpha=1,zorder=1)

                    #Note: Add highlight for specific time period (optional)
                    hli1 = dt.datetime(2002,3,18,14,50,15)
                    hlf1 = dt.datetime(2002,3,18,14,50,55)
                    hli2 = dt.datetime(2002,3,18,14,51,10)
                    hlf2 = dt.datetime(2002,3,18,14,51,40)
                    ax[i,0].axvspan(hli1,hlf1,facecolor='lightgray',edgecolor='None',alpha=0.8, zorder=0)
                    ax[i,0].axvspan(hli2,hlf2,facecolor='lightgray',edgecolor='None',alpha=0.8, zorder=0)
                    
                #Set x-axis properties
                dtime_ticks = []
                dtime_ticklabels = []
                if (self.dtf-self.dti).seconds > 1200:
                    major_tick_interval = 600
                else:
                    major_tick_interval = 60 
                dtime_range = (self.dtf - self.dti).seconds/major_tick_interval
                m = 0
                while m <= dtime_range:
                    tick_dt = self.dti + m*dt.timedelta(seconds=major_tick_interval)
                    dtime_ticks.append(tick_dt)
                    tick_time_split = str(tick_dt).split(' ')[1].split(':')
                    tick_dt_str = ':'.join(
                        (tick_time_split[0], tick_time_split[1]))
                    dtime_ticklabels.append(tick_dt_str)
                    m = m+1

                if xax_view == 0:
                    if i < no_subplots-1:
                        ax[i,0].minorticks_on()
                        ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                        ax[i,0].tick_params(
                            axis='x', which='both', bottom=True, top=False, labelbottom=False)
                        ax[i,0].tick_params(
                            axis='x', which='major', direction='inout', width=1, length=15)
                        ax[i,0].axes.get_xaxis().set_ticks(
                            dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize)
                        ax[i,0].tick_params(
                            axis='x', which='minor', direction='inout', width=1, length=5)
                        ax[i,0].set_xlim(self.dti, self.dtf)
                    else:
                        ax[i,0].minorticks_on()
                        ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                        ax[i,0].tick_params(
                            axis='x', which='both', bottom=True, top=False, labelbottom=True)
                        ax[i,0].tick_params(
                            axis='x', which='major', direction='inout', width=1, length=15)
                        ax[i,0].axes.get_xaxis().set_ticks(
                            dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize)
                        ax[i,0].tick_params(
                            axis='x', which='minor', direction='inout', width=1, length=5)
                        ax[i,0].set_xlim(self.dti, self.dtf)
                        dt_str = self.dti.strftime(
                            '%Y-%m-%d')
                        ax[i,0].set_xlabel('{} (UT)'.format(
                            dt_str), fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif xax_view == 1:
                    ax[i,0].minorticks_on()
                    ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                    ax[i,0].tick_params(
                        axis='x', which='both', bottom=True, top=False, labelbottom=False)
                    ax[i,0].tick_params(
                        axis='x', which='major', direction='inout', width=1, length=15)
                    ax[i,0].axes.get_xaxis().set_ticks(
                        dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize) #must set to ensure ticks are in correct position
                    ax[i,0].tick_params(
                        axis='x', which='minor', direction='inout', width=1, length=5)
                    ax[i,0].set_xlim(self.dti, self.dtf)
                elif xax_view == 2:
                    ax[i,0].minorticks_on()
                    ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                    ax[i,0].tick_params(
                        axis='x', which='both', bottom=True, top=False, labelbottom=True)
                    ax[i,0].tick_params(
                        axis='x', which='major', direction='inout', width=1, length=15)
                    ax[i,0].axes.get_xaxis().set_ticks(
                        dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize)
                    ax[i,0].tick_params(
                        axis='x', which='minor', direction='inout', width=1, length=5)
                    ax[i,0].set_xlim(self.dti, self.dtf)
                    dt_str = self.dti.strftime(
                        '%Y-%m-%d')  # set date string
                    ax[i,0].set_xlabel('{} (UT)'.format(
                        dt_str), fontsize=axlabel_fontsize, fontfamily='sans-serif')
                    
                #Set y-axis properties
                ax[i,0].yaxis.set_minor_locator(AutoMinorLocator(5))
                ax[i,0].tick_params(axis='y', which='both',
                                    left=True, right=False)
                ax[i,0].tick_params(axis='y', which='major',
                                    direction='inout', width=1, length=15, labelsize=axtick_fontsize)
                ax[i,0].tick_params(axis='y', which='minor',
                                    direction='inout', width=1, length=5)
                
                    
                #Set y-axis limits and other properties.
                #Set zero for minimum limit.
                #Maximum limit depends on maximum values involved.
                ymin = 0 
                ymax = 1 #Set extreme values that will be satisfied in first iteration of loop below.
                
                #Loop over each summed flux parameter. If maximum value of array is above current limit, update limit.
                for key in PanelDict['LineSum'][self.sc].keys():
                    if np.max(np.min(PanelDict['LineSum'][self.sc][key][:,1])) > ymax:
                        ymax = np.max(PanelDict['LineSum'][self.sc][key][:,1])
                        ymax = 1.2*ymax
                ax[i,0].set_ylim([ymin, ymax])
                ax[i,0].set_ylabel('Flux >{} keV'.format(self.energy_min) + '\n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$', 
                                   fontsize=axlabel_fontsize, fontfamily='sans-serif')
                ax[i,0].ticklabel_format(axis='y',style='scientific',scilimits=(3,3))
                    

                #Set in-plot text labels
                #Panel letters
                ax[i,0].text(0.04, 0.85, '({})'.format(ALC[i+self.sp_ind0]), fontsize=subtitle_fontsize, fontfamily='sans-serif',
                             horizontalalignment='center', verticalalignment='center',
                             transform=ax[i,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                #Spacecraft indicators
                ax[i,0].text(0.96, 0.85, self.sc, fontsize=subtitle_fontsize, fontfamily='sans-serif',
                             horizontalalignment='center', verticalalignment='center',
                             transform=ax[i,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                #Set legend.
                ax[i,0].legend(loc='lower right', prop={'size': legend_fontsize,'family': 'sans-serif'})
                
                #Set axes colors depending on light or dark modes.
                ax[i,0].spines['bottom'].set_color(ctick)
                ax[i,0].spines['top'].set_color(ctick)
                ax[i,0].spines['left'].set_color(ctick)
                ax[i,0].spines['right'].set_color(ctick)
                ax[i,0].tick_params(axis='both',which='both',colors=ctick)
                ax[i,0].xaxis.label.set_color(ctick) 
                ax[i,0].yaxis.label.set_color(ctick) 
                
                #Save plot to CSV file
                if data_record == True:
                    #Set CSV file name.
                    start_dtstr = self.dti.strftime('%Y-%m-%dT%H-%M-%SZ')
                    stop_dtstr = self.dtf.strftime('%Y-%m-%dT%H-%M-%SZ')
                    csv_filename = '_'.join(['LineSum',self.sc,start_dtstr,stop_dtstr]) + '.csv'
                    
                    #Open CSV file.
                    with open(csv_filename,'w', newline = '') as csvfile:
                        csvwriter = csv.writer(csvfile)
                    
                        #All parameters for the same Cluster spacecraft share a datetime list.
                        #Use the first one
                        para0 = list(PanelDict['LineSum'][self.sc].keys())[0]
                        dt_arr = PanelDict['LineSum'][self.sc][para0][:,0]
                        
                        #Each datetime represents a row in the CSV file.
                        for i in np.arange(len(dt_arr)):
                            
                            #Build a row containing all data
                            row_data = [dt_arr[i]]
                            for key in PanelDict['LineSum'][self.sc].keys():
                                data = PanelDict['LineSum'][self.sc][key][i,1]
                                row_data.append(data)
                            
                            csvwriter.writerow(row_data)
                    csvfile.close()
                    print('LineSum, {} data recorded in csv'.format(self.sc))

        #Save figure
        if save == True:
            targetfolder = 'ClusterPanels/{}/{}/{}/'.format(self.dti.year,self.dti.month,self.dti.day)
            os.makedirs(targetfolder,exist_ok=True)
            start_str = dt.datetime.strftime(self.dti, '%Y%m%d%H%M%S')
            stop_str = dt.datetime.strftime(self.dtf, '%Y%m%d%H%M%S')
            savename = targetfolder + '{}_{}_{}'.format(start_str,stop_str,self.sc)
    
            for para in paralist:
                savename = savename + '_{}'.format(para)
    
            savename = savename + '.{}'.format(saveformat)
            plt.subplots_adjust(left=0.12, bottom=0.05, right=0.90, top=0.97,
                                wspace=0.0, hspace=0.08)
    
            plt.savefig(savename, dpi=200, format=saveformat, transparent=True)
            plt.close()


    def csv_record(self,para,sc,comb_data):
        #Method for recording data to CSV file (Not for LineSum).
        start_dtstr = self.dti.strftime('%Y-%m-%dT%H-%M-%SZ')
        stop_dtstr = self.dtf.strftime('%Y-%m-%dT%H-%M-%SZ')
        csv_filename = '_'.join([para,sc,start_dtstr,stop_dtstr]) + '.csv'
        with open(csv_filename,'w', newline = '') as csvfile:
            csvwriter = csv.writer(csvfile)
            
            for i in np.arange(len(comb_data[:,0])):
                csvwriter.writerow(comb_data[i,:])
                            
        csvfile.close()

"""
Operating Script (Figure S3)
"""
"""
dti = dt.datetime(2002,3,18,14,49,0)
dtf = dt.datetime(2002,3,18,14,54,0)
paralist = ['H+_fluxU_Sun','H+_fluxU_Dawn','H+_fluxU_Dusk','H+_fluxU_Antisun']
sc = 'C4'
energy_min = 10 #in keV
FLS = FluxLineSum(dti,dtf,paralist,sc,energy_min,new_dl=True,just_sum=False,autoplot=True,label_start='A',
             xax_view=0,mode='light',save=True,saveformat='png',data_record=True)
"""