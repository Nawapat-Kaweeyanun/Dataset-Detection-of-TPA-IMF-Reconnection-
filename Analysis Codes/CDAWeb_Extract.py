"""
Edition Date: 2025-October-20
@author: Nawapat Kaweeyanun
"""
"""
Objective: Download and plot OMNI data from NASA's Coordinate Data Analysis (CDAS/CDAWeb)

Prerequisit: cdasws module

Note: Example script for application provided at the end.
"""


"""
Define class for CDA data manager
"""

from cdasws import CdasWs
from cdasws.datarepresentation import DataRepresentation
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from string import ascii_uppercase as ALC
from matplotlib.ticker import AutoMinorLocator
from spacepy import coordinates as coord
import csv
import os


class CDAWeb_Manager(object):

    def __init__(self, start_dt, end_dt, paralist, observatoryGroup, instrumentType=None, res='1min', grid='HRO2', coordinates ='GSE'):        
        self.start_dt = start_dt #Period start datetime
        self.end_dt = end_dt #Period end datetime
        self.OG = observatoryGroup #Mission to obtain (e.g., OMNI)
        self.IT = instrumentType #Instrument to analyze (e.g., magnetometer)
        self.paralist = paralist #List of parameters (e.g., magnetic field/velocity components)
        self.res = res #Data resolution (1 hr, 5 min, or 1 min)
        self.grid = grid #Grid (HRO (old, not recommended) or HRO2 (new, recommended))
        self.coord = coord #Coordinate (GSE or GSM)

        #Call CdasWs module
        self.cdas = CdasWs()
        
        #Call OMNI acquisition method after checking inputs are correct
        if 'OMNI' in self.OG:
        
            if self.res not in ['1hr', '5min', '1min']:
                raise Exception('Resolution incorrect')
        
            if self.grid not in ['HRO2', 'HRO']:
                raise Exception('Grid type incorrect')
        
            self.OMNI_Acquire(self.res, self.grid, self.coord)

    def OMNI_Acquire(self, res, grid):

        #Get OMNI data from CdasWs module object
        datasets = self.cdas.get_datasets(observatoryGroup=self.OG,
                                          instrumentType=self.IT)

        #Make sure the input coordinates are correct
        if coord not in ['GSE', 'GSM']:
            raise Exception('Coordinate type incorrect')

        #Filter data for matching time resolution and grid
        for n in np.arange(len(datasets)):
            ID = datasets[n]['Id']
            if res.upper() in ID and grid.upper() in ID:
                MatchID = ID
                break
            else:
                MatchID = None
                continue

        if MatchID == None:
            raise Exception(
                'Dataset not found for OMNI at res {} and grid {}'.format(res, grid))

        #Convert input parameter list into OMNI-appropriate parameter IDs
        self.varlist = []
        for para in self.paralist:
            if res == '1hr':
                if para in ['F','V','N','T','E','Pressure','Beta','THETA-V','PHI-V']:
                    self.varlist.append('{}_1800').format(para) #No upper and case sensitive
                elif para in ['Bx']:
                    self.varlist.append('BX_GSE1800') #Bx only in GSE
                elif para in ['By','Bz']:                
                    self.varlist.append('{}_{}1800'.format(para.upper(),self.coord))
            else:
                if para in ['F','Vx','Vy','Vz','x','y','z','T','E','Pressure','Beta']:
                    self.varlist.append('{}'.format(para))
                elif para in ['Bx']:
                    self.varlist.append('BX_GSE')
                elif para in ['By','Bz']:
                    self.varlist.append('{}_{}'.format(para.upper(),self.coord))
                elif para in ['V']:
                    self.varlist.append('flow_speed')
                elif para in ['N']:
                    self.varlist.append('proton_density')
        
        #Obtain ordered data as xarrays
        self.data = self.cdas.get_data(MatchID, self.varlist, self.start_dt, self.end_dt,
                                       dataRepresentation=DataRepresentation.XARRAY)[1]
        
    def OMNI_CA(self):
        #Calculate clock angle if required.
        #Clock angle: like a clock if look from the Sun i.e., 90 along YGSE > 0 axis
        
        #Check that clock angle is specified in parameter list
        if 'CA' not in self.paralist:
            raise Exception('Clock Angle not specified in parameter list.')

        #Search varlist for By and Bz (will only have one version depending on res and coord)
        By_varname = [var for var in self.varlist if 'BY' in var][0] #call 0th element to get string
        Bz_varname = [var for var in self.varlist if 'BZ' in var][0]
        
        #Obtain data array
        By_arr = self.data[By_varname].values
        Bz_arr = self.data[Bz_varname].values
        
        #Calculate clock angles
        #Note: No need to convert ang = 2pi-ang because arctan2 going anticlockwise when viewed toward the Sun = clock angle going clockwise when viewed from the Sun.
        CA_arr = np.arctan2(By_arr,Bz_arr) #Note: it's atan2(y,x). Bz is x in this case because reference axis is vertical.
        CA_arr = (CA_arr)*180/np.pi #Convert radian to degrees
        
        #Include CA_arr into self and adds to paralist and varlist
        self.CA = CA_arr
        self.varlist.append('CA')


    def OMNI_Plot(self, xax_view=0, label_start='A', veccomb=False, saveformat='png',
                  MP_HL=False, Walen_HL=False, mode='light', data_record=False):
        #Plot OMNI data from xarrays
        #If combine vector components, create new para list of new length (placeholder) where
        #components are combined (but only if there are all three, otherwise plotting choice gets too hard later)
        
        paratemp = self.paralist.copy()
        varlist = self.varlist.copy()

        #Readjust variable list in case vector components are combined.
        if veccomb == True:
            if 'Bx' in paratemp:
                paratemp.remove('Bx')
                #Have to adjust varlist too in case non Bvec/Vvec data are put in so plots can work
                varlist.remove([var for var in self.varlist if 'BX' in var][0])
                try:  #Remove the other two if possible
                    paratemp.remove('By')
                    varlist.remove([var for var in self.varlist if 'BY' in var][0])
                    paratemp.remove('Bz')
                    varlist.remove([var for var in self.varlist if 'BZ' in var][0])
                except:
                    raise Exception(
                        'Not all three components present for combination')
                paratemp.append('Bvec')
            if 'Vx' in paratemp:
                paratemp.remove('Vx')
                varlist.remove([var for var in self.varlist if 'Vx' in var][0])
                try:  # remove the other two if possible
                    paratemp.remove('Vy')
                    varlist.remove([var for var in self.varlist if 'Vy' in var][0])
                    paratemp.remove('Vz')
                    varlist.remove([var for var in self.varlist if 'Vz' in var][0])
                except:
                    raise Exception(
                        'Not all three components present for combination')
                paratemp.append('Vvec')

        
        #Make figure with appropriate number of subplots
        fig, ax = plt.subplots(len(paratemp), 1, figsize=(20, 5*len(paratemp)),squeeze=False)
        subtitle_fontsize = 32
        axlabel_fontsize = 24
        axtick_fontsize = 18
                
        #Define figure properties for light and dark mode
        if mode == 'light':
            ctick = 'k' #Black ticks/labels for scalars
            cline = 'k' 
            cxline = 'r' #Set color for x-y-z vectors
            cyline = 'k'
            czline = 'b'
            tchoice = False #Plot with white background
        elif mode == 'dark':
            ctick = 'w' #White ticks/labels for scalars
            cline = 'w'
            cxline = 'r' #Set color for x-y-z vectors
            cyline = 'w'
            czline = 'cyan'
            tchoice = True #Plot with transparent background
        else:
            raise Exception('Incorrect Colour Mode')

        #Set limits and labels for different variables
        def lim_unit_setter(var):
        
            if var in ['F1800', 'F']:
                ymin = 0
                ymax = 60
                unit = '|B| (nT)'
            elif var in ['BX_GSE1800', 'BX_GSE', 'BY_GSE1800', 'BY_GSE', 'BZ_GSE1800', 'BZ_GSE',
                         'BX_GSM1800', 'BX_GSM', 'BY_GSM1800', 'BY_GSM', 'BZ_GSM1800', 'BZ_GSM']:
                ymin = -25
                ymax = 25
                unit = 'B{} (nT)'.format(var[1].lower())
            elif var in ['Vx', 'Vy', 'Vz']:
                ymin = -800
                ymax = 200
                unit = 'V{} (km/s)'.format(var[1].lower())
            elif var == 'CA':
                ymin = -180
                ymax = 180
                unit = 'Clock Angle (deg)'
            else:
                ymin = np.amin(self.data[var]) - 0.2 * \
                    abs(np.amin(self.data[var]))  # +/- 20%
                ymax = np.amax(self.data[var]) + 0.2 * \
                    abs(np.amax(self.data[var]))
                unit = ''
            return ymin, ymax, unit

        if veccomb == False: #Vector components separated into multiple panels.

            for n in np.arange(len(varlist)):
                
                #Special clause for clock angle
                if varlist[n] == 'CA':
                    ax[n,0].plot(self.data['Epoch'].values, self.CA, color=cline,lw=2)
                
                else:
                    ax[n,0].plot(self.data['Epoch'].values,
                               self.data[varlist[n]].values, color=cline,lw=2)

                #If preferred, highlight area under curve for Bz
                #if self.varlist[n] == 'BZ_GSE':
                #    ax[n].fill_between(self.data['Epoch'].values,self.data[self.varlist[n]].values,np.zeros(len(self.data['Epoch'])))
                
                #Produce datetime string for x-axis labelling
                dt64 = self.data['Epoch'].values[0]
                unix_epoch = np.datetime64(0, 's')
                one_second = np.timedelta64(1, 's')
                seconds_since_epoch = (dt64-unix_epoch)/one_second
                dt_start = dt.datetime.utcfromtimestamp(seconds_since_epoch)
                dtstr = dt_start.strftime('%Y-%m-%d')

                #Produce dt.datetime for end datetime (for magnetopause highlight below)
                dt_end = dt.datetime.utcfromtimestamp(
                    (self.data['Epoch'].values[-1]-unix_epoch)/one_second)

                #Highlight magnetopause crossing time
                if MP_HL == True:
                    MP_start = dt.datetime(2002, 3, 18, 14, 54, 52)
                    MP_end = dt.datetime(2002, 3, 18, 15, 3, 52)

                    #Check that MP crossing time is within data datetime range
                    if MP_start < dt_start or MP_end > dt_end:
                        print(
                            'Magnetopause crossing interval not in plot range. No highlight.')
                    else:
                        ax[n,0].axvspan(
                                MP_start, MP_end, facecolor='lightgray', edgecolor='None', alpha=0.5, zorder=0)
                
                #Highlight Walen test interval
                if Walen_HL == True:
                    Walen_start = dt.datetime(2002, 3, 18, 14, 57, 25)
                    Walen_end = dt.datetime(2002, 3, 18, 15, 3, 6)

                    if Walen_start < dt_start or Walen_end > dt_end:
                        print('Walen Test interval not in plot range. No highlight.')
                    else:
                        ax[n,0].axvline(Walen_start, color='purple',
                                        ls='--', lw=1.5, zorder=1)
                        ax[n,0].axvline(Walen_end, color='purple',
                                        ls='--', lw=1.5, zorder=1)
                
                #x-axis label settings
                dtime_ticks = []
                dtime_ticklabels = []
                start_dtobj = dt.datetime.strptime(self.start_dt,'%Y-%m-%dT%H:%M:%SZ')
                end_dtobj = dt.datetime.strptime(self.end_dt,'%Y-%m-%dT%H:%M:%SZ')
                
                #10 minutes major tick intervals for time period more than 20 minutes, 1 minute intervals for less.
                if (end_dtobj-start_dtobj).seconds > 1200:
                    major_tick_interval = 600 #seconds
                else:
                    major_tick_interval = 60 #seconds
                dtime_range = (end_dtobj - start_dtobj).seconds/major_tick_interval
                m = 0
                while m <= dtime_range:
                    tick_dt = start_dtobj + m*dt.timedelta(seconds=major_tick_interval)
                    dtime_ticks.append(tick_dt)
                    tick_time_split = str(tick_dt).split(' ')[1].split(':')
                    tick_dt_str = ':'.join(
                        (tick_time_split[0], tick_time_split[1]))
                    dtime_ticklabels.append(tick_dt_str)
                    m = m+1
                
                #set axis settings
                ax[n,0].set_xlim([self.data['Epoch'].values[0],
                               self.data['Epoch'].values[-1]])
                if xax_view <= 0:  #Only the bottom-most subplot has x-axis label
                    if n < len(paratemp)-1:
                        ax[n,0].minorticks_on()
                        ax[n,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                        ax[n,0].tick_params(
                            axis='x', which='both', bottom=True, top=False, labelbottom=False)
                        ax[n,0].tick_params(
                            axis='x', which='major', direction='inout', width=1, length=15)
                        ax[n,0].axes.get_xaxis().set_ticks(
                            dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize)
                        ax[n,0].tick_params(
                            axis='x', which='minor', direction='inout', width=1, length=5)
                        ax[n,0].set_xlim(self.data['Epoch'].values[0],
                                       self.data['Epoch'].values[-1])
                    else:
                        ax[n,0].minorticks_on()
                        ax[n,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                        ax[n,0].tick_params(
                            axis='x', which='both', bottom=True, top=False, labelbottom=True)
                        ax[n,0].tick_params(
                            axis='x', which='major', direction='inout', width=1, length=15)
                        ax[n,0].axes.get_xaxis().set_ticks(
                            dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize)
                        ax[n,0].tick_params(
                            axis='x', which='minor', direction='inout', width=1, length=5)
                        ax[n,0].set_xlim(self.data['Epoch'].values[0],
                                       self.data['Epoch'].values[-1])
                        ax[n,0].set_xlabel(dtstr, fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif xax_view == 1:  #No subplot has x-axis label
                    ax[n,0].minorticks_on()
                    ax[n,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                    ax[n,0].tick_params(axis='x', which='both',
                                      bottom=True, top=False, labelbottom=False)
                    ax[n,0].tick_params(axis='x', which='major',
                                      direction='inout', width=1, length=15)
                    ax[n,0].axes.get_xaxis().set_ticks(
                        dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize)
                    ax[n,0].tick_params(axis='x', which='minor',
                                      direction='inout', width=1, length=5)
                    ax[n,0].set_xlim(self.data['Epoch'].values[0],
                                   self.data['Epoch'].values[-1])
                elif xax_view >= 2:  #All subplots have x-axis label
                    ax[n,0].minorticks_on()
                    ax[n,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                    ax[n,0].tick_params(axis='x', which='both',
                                      bottom=True, top=False, labelbottom=False)
                    ax[n,0].tick_params(axis='x', which='major',
                                      direction='inout', width=1, length=15)
                    ax[n,0].axes.get_xaxis().set_ticks(
                        dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize)
                    ax[n,0].tick_params(axis='x', which='minor',
                                      direction='inout', width=1, length=5)
                    ax[n,0].set_xlim(self.data['Epoch'].values[0],
                                   self.data['Epoch'].values[-1])
                    ax[n,0].xaxis.set_major_formatter(
                        mdates.DateFormatter("%H:%M"))
                    ax[n,0].set_xlabel(dtstr, fontsize=axlabel_fontsize, fontfamily='sans-serif')

                #y-axis label settings
                ymin, ymax, unit = lim_unit_setter(self.varlist[n])
                ax[n,0].yaxis.set_minor_locator(AutoMinorLocator(5))
                ax[n,0].tick_params(axis='y', which='major',
                                  direction='inout', width=1, length=15)
                ax[n,0].tick_params(axis='y', which='minor',
                                  direction='inout', width=1, length=5)
                ax[n,0].set_ylim(ymin, ymax)
                ax[n,0].set_ylabel(unit, fontsize=axlabel_fontsize, fontfamily='sans-serif')
                ax[n,0].tick_params(axis='both', which='major',
                                  labelsize=axtick_fontsize)  
                #Add zero line if needed
                if ymin < 0 and ymax > 0:
                    ax[n,0].axhline(y=0,color=ctick,lw=1,alpha=1,zorder=1)
                

                #Set in-plot labels
                sp_ind0 = ALC.find(label_start.upper())
                ax[n,0].text(0.04, 0.85, '({})'.format(ALC[n+sp_ind0]), fontsize=subtitle_fontsize, fontfamily='sans-serif',
                           horizontalalignment='center', verticalalignment='center',
                           transform=ax[n,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                ax[n,0].text(0.95, 0.85, 'OMNI', fontsize=subtitle_fontsize, fontfamily='sans-serif',
                           horizontalalignment='center', verticalalignment='center',
                           transform=ax[n,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                
                #Set colors
                ax[n,0].spines['bottom'].set_color(ctick) #Set axes color
                ax[n,0].spines['top'].set_color(ctick)
                ax[n,0].spines['left'].set_color(ctick)
                ax[n,0].spines['right'].set_color(ctick)
                ax[n,0].tick_params(axis='both',which='both',colors=ctick) #Set tick color
                ax[n,0].xaxis.label.set_color(ctick) #Set x-axis ticklabel color
                ax[n,0].yaxis.label.set_color(ctick) #Set y-axis ticklabel color  

                #If needed, record plot data in CSV file
                if data_record == True:
                    self.csv_record(varlist[n])

        elif veccomb == True: #Vector components combined

            #Set label of first subplot (for publication)
            sp_ind0 = ALC.find(label_start.upper())
            

            for n in np.arange(len(paratemp)):
                para = paratemp[n]

                if para == 'Bvec':
                    
                    #Identify correct Bx, By, Bz variable names.
                    #Only one variant of Bx, By, Bz will exist in the variable list.
                    Bx_varname = [var for var in self.varlist if 'BX' in var][0]
                    By_varname = [var for var in self.varlist if 'BY' in var][0]
                    Bz_varname = [var for var in self.varlist if 'BZ' in var][0]
                    
                    #Plot three vector components in the same subplot.
                    #Otherwise, settings same as above.
                    ax[n,0].plot(self.data['Epoch'].values, self.data[Bx_varname].values,
                               label='X', color=cxline, ls='-', lw=2, alpha=1, zorder=2)
                    ax[n,0].plot(self.data['Epoch'].values, self.data[By_varname].values,
                               label='Y', color=cyline, ls='-', lw=2, alpha=1, zorder=2)
                    ax[n,0].plot(self.data['Epoch'].values, self.data[Bz_varname].values,
                               label='Z', color=czline, ls='-', lw=2, alpha=1, zorder=2)
                    
                    #y-axis settings
                    ymin, ymax, _ = lim_unit_setter(Bx_varname)
                    ax[n,0].minorticks_on() 
                    ax[n,0].set_ylim(ymin, ymax)
                    ax[n,0].yaxis.set_minor_locator(AutoMinorLocator(5))
                    ax[n,0].tick_params(
                        axis='y', which='major', direction='inout', width=1, length=15)
                    ax[n,0].tick_params(
                        axis='y', which='minor', direction='inout', width=1, length=5)
                    ax[n,0].set_ylabel(
                        '$\mathrm{B_{IMF}}$ (nT)', fontsize=axlabel_fontsize, fontfamily='sans-serif',)
                    ax[n,0].tick_params(
                        axis='both', which='major', labelsize=axtick_fontsize)

                    if ymin < 0 and ymax > 0:
                        ax[n,0].axhline(y=0,color=ctick,lw=1,alpha=1,zorder=1)
                    ax[n,0].text(0.04, 0.85, '({})'.format(ALC[n+sp_ind0]), fontsize=subtitle_fontsize, fontfamily='sans-serif',
                               horizontalalignment='center', verticalalignment='center',
                               transform=ax[n,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                    ax[n,0].text(0.95, 0.85, 'OMNI', fontsize=subtitle_fontsize, fontfamily='sans-serif',
                               horizontalalignment='center', verticalalignment='center',
                               transform=ax[n,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                    
                    ax[n,0].spines['bottom'].set_color(ctick) 
                    ax[n,0].spines['top'].set_color(ctick)
                    ax[n,0].spines['left'].set_color(ctick)
                    ax[n,0].spines['right'].set_color(ctick)
                    ax[n,0].tick_params(axis='both',which='both',colors=ctick)
                    ax[n,0].xaxis.label.set_color(ctick) 
                    ax[n,0].yaxis.label.set_color(ctick) 

                elif para == 'Vvec': #For solar wind velocity, which has specific variable names.
                    ax[n,0].plot(self.data['Epoch'].values, self.data['Vx'].values,
                               label='Vx', color=cxline, ls='-', lw=2, alpha=1, zorder=2)
                    ax[n,0].plot(self.data['Epoch'].values, self.data['Vy'].values,
                               label='Vy', color=cyline, ls='-', lw=2, alpha=1, zorder=2)
                    ax[n,0].plot(self.data['Epoch'].values, self.data['Vz'].values,
                               label='Vz', color=czline, ls='-', lw=2, alpha=1, zorder=2)
                    ymin, ymax, _ = lim_unit_setter('Vx')
                    ax[n,0].minorticks_on()
                    ax[n,0].set_ylim(ymin, ymax)
                    ax[n,0].yaxis.set_minor_locator(AutoMinorLocator(5))
                    ax[n,0].tick_params(axis='y', which='major',
                                      direction='inout', width=1, length=15)
                    ax[n,0].tick_params(axis='y', which='minor',
                                      direction='inout', width=1, length=5)
                    ax[n,0].set_ylabel(
                        '$\mathrm{V_{SW,GSE}}$ (nT)', fontsize=axlabel_fontsize, fontfamily='sans-serif')
                    ax[n,0].tick_params(axis='both', which='major',
                                      labelsize=axtick_fontsize)
                    if ymin < 0 and ymax > 0:
                        ax[n,0].axhline(y=0,color=ctick,lw=1,alpha=1,zorder=1)
                    ax[n,0].text(0.04, 0.85, '({})'.format(ALC[n+sp_ind0]), fontsize=subtitle_fontsize, fontfamily='sans-serif',
                               horizontalalignment='center', verticalalignment='center',
                               transform=ax[n,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                    ax[n,0].text(0.95, 0.85, 'OMNI', fontsize=subtitle_fontsize,
                               horizontalalignment='center', verticalalignment='center', fontfamily='sans-serif',
                               transform=ax[n,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))

                    ax[n,0].spines['bottom'].set_color(ctick)
                    ax[n,0].spines['top'].set_color(ctick)
                    ax[n,0].spines['left'].set_color(ctick)
                    ax[n,0].spines['right'].set_color(ctick)
                    ax[n,0].tick_params(axis='both',which='both',colors=ctick)
                    ax[n,0].xaxis.label.set_color(ctick) 
                    ax[n,0].yaxis.label.set_color(ctick)  
                
                else:  #Other parameters (scalar)
                
                    #Special clause for clock angle plotting
                    if self.varlist[n] == 'CA': 
                        ax[n,0].plot(self.data['Epoch'].values, self.CA, color=cline, lw=2,  alpha=1, zorder=2)
                    else:
                        ax[n,0].plot(self.data['Epoch'].values, self.data[self.varlist[n]].values, 
                                     color=cline, lw=2,  alpha=1, zorder=2)
                    ymin, ymax, unit = lim_unit_setter(self.varlist[n])
                    ax[n,0].minorticks_on()
                    ax[n,0].set_ylim(ymin, ymax)
                    ax[n,0].yaxis.set_minor_locator(AutoMinorLocator(5))
                    ax[n,0].tick_params(axis='y', which='major',
                                      direction='inout', width=1, length=15)
                    ax[n,0].tick_params(axis='y', which='minor',
                                      direction='inout', width=1, length=5)
                    ax[n,0].set_ylabel(unit, fontsize=axlabel_fontsize, fontfamily='sans-serif')
                    ax[n,0].tick_params(axis='both', which='major',
                                      labelsize=axtick_fontsize)
                    ax[n,0].text(0.04, 0.85, '({})'.format(ALC[n+sp_ind0]), fontsize=subtitle_fontsize, fontfamily='sans-serif',
                               horizontalalignment='center', verticalalignment='center',
                               transform=ax[n,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                    ax[n,0].text(0.95, 0.85, 'OMNI', fontsize=subtitle_fontsize, fontfamily='sans-serif',
                               horizontalalignment='center', verticalalignment='center',
                               transform=ax[n,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                    

                    ax[n,0].spines['bottom'].set_color(ctick)
                    ax[n,0].spines['top'].set_color(ctick)
                    ax[n,0].spines['left'].set_color(ctick)
                    ax[n,0].spines['right'].set_color(ctick)
                    ax[n,0].tick_params(axis='both',which='both',colors=ctick) 
                    ax[n,0].xaxis.label.set_color(ctick) 
                    ax[n,0].yaxis.label.set_color(ctick)

                #Produce datetime string for x labelling (use start of range)
                dt64 = self.data['Epoch'].values[0]
                unix_epoch = np.datetime64(0, 's')
                one_second = np.timedelta64(1, 's')
                seconds_since_epoch = (dt64-unix_epoch)/one_second
                dt_start = dt.datetime.utcfromtimestamp(seconds_since_epoch)
                dtstr = dt_start.strftime('%Y-%m-%d')

                #Produce dt.datetime for end datetime (for highlights)
                dt_end = dt.datetime.utcfromtimestamp(
                    (self.data['Epoch'].values[-1]-unix_epoch)/one_second)

                #Highlight MP crossing times
                if MP_HL == True:
                    MP_start = dt.datetime(2002, 3, 18, 14, 54, 52)
                    MP_end = dt.datetime(2002, 3, 18, 15, 3, 52)
                    

                    if MP_start < dt_start or MP_end > dt_end:
                        print(
                            'Magnetopause crossing interval not in plot range. No highlight.')
                    else:
                        ax[n,0].axvspan(
                                MP_start, MP_end, facecolor='lightgray', edgecolor='None', alpha=0.5, zorder=0)
                
                #Highlight Walen test interval
                if Walen_HL == True:
                    Walen_start = dt.datetime(2002, 3, 18, 14, 57, 25)
                    Walen_end = dt.datetime(2002, 3, 18, 15, 3, 6)

                    if Walen_start < dt_start or Walen_end > dt_end:
                        print('Walen Test interval not in plot range. No highlight.')
                    else:
                        ax[n,0].axvline(Walen_start, color='purple',
                                        ls='--', lw=1.5, zorder=1)
                        ax[n,0].axvline(Walen_end, color='purple',
                                        ls='--', lw=1.5, zorder=1)


                #Set legend
                ax[n,0].legend(loc='lower right', prop={'size': 14, 'family':'sans-serif'})

                #x-axis settings
                ax[n,0].set_xlim([self.data['Epoch'].values[0],
                               self.data['Epoch'].values[-1]])
                
                dtime_ticks = []
                dtime_ticklabels = []          
                start_dtobj = dt.datetime.strptime(self.start_dt,'%Y-%m-%dT%H:%M:%SZ')
                end_dtobj = dt.datetime.strptime(self.end_dt,'%Y-%m-%dT%H:%M:%SZ')            
                if (end_dtobj-start_dtobj).seconds > 1200:
                    major_tick_interval = 600 
                else:
                    major_tick_interval = 60 
                dtime_range = (end_dtobj - start_dtobj).seconds/major_tick_interval
                m = 0
                while m <= dtime_range:
                    tick_dt = start_dtobj + m*dt.timedelta(seconds=major_tick_interval)
                    dtime_ticks.append(tick_dt)
                    tick_time_split = str(tick_dt).split(' ')[1].split(':')
                    tick_dt_str = ':'.join(
                        (tick_time_split[0], tick_time_split[1]))
                    dtime_ticklabels.append(tick_dt_str)
                    m = m+1
                
                if xax_view <= 0:
                    if n < len(paratemp)-1:
                        ax[n,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                        ax[n,0].tick_params(
                            axis='x', which='both', bottom=True, top=False, labelbottom=False)
                        ax[n,0].tick_params(
                            axis='x', which='major', direction='inout', width=1, length=15)
                        ax[n,0].axes.get_xaxis().set_ticks(
                            dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize)
                        ax[n,0].tick_params(
                            axis='x', which='minor', direction='inout', width=1, length=5)
                        ax[n,0].set_xlim(self.data['Epoch'].values[0],
                                       self.data['Epoch'].values[-1])
                    else:
                        ax[n,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                        ax[n,0].tick_params(
                            axis='x', which='both', bottom=True, top=False, labelbottom=True)
                        ax[n,0].tick_params(
                            axis='x', which='major', direction='inout', width=1, length=15)
                        ax[n,0].axes.get_xaxis().set_ticks(
                            dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize)
                        ax[n,0].tick_params(
                            axis='x', which='minor', direction='inout', width=1, length=5)
                        ax[n,0].set_xlim(self.data['Epoch'].values[0],
                                       self.data['Epoch'].values[-1])
                        ax[n,0].xaxis.set_major_formatter(
                            mdates.DateFormatter("%H:%M"))
                        ax[n,0].set_xlabel(dtstr, fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif xax_view == 1:
                    ax[n,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                    ax[n,0].tick_params(axis='x', which='both',
                                      bottom=True, top=False, labelbottom=False)
                    ax[n,0].tick_params(axis='x', which='major',
                                      direction='inout', width=1, length=15)
                    ax[n,0].axes.get_xaxis().set_ticks(
                        dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize)
                    ax[n,0].tick_params(axis='x', which='minor',
                                      direction='inout', width=1, length=5)
                    ax[n,0].set_xlim(self.data['Epoch'].values[0],
                                   self.data['Epoch'].values[-1])
                elif xax_view >= 2:
                    ax[n,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                    ax[n,0].tick_params(axis='x', which='both',
                                      bottom=True, top=False, labelbottom=True)
                    ax[n,0].tick_params(axis='x', which='major',
                                      direction='inout', width=1, length=15)
                    ax[n,0].axes.get_xaxis().set_ticks(
                        dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize)
                    ax[n,0].tick_params(axis='x', which='minor',
                                      direction='inout', width=1, length=5)
                    ax[n,0].set_xlim(self.data['Epoch'].values[0],
                                   self.data['Epoch'].values[-1])
                    ax[n,0].xaxis.set_major_formatter(
                        mdates.DateFormatter("%H:%M"))
                    ax[n,0].set_xlabel(dtstr, fontsize=axlabel_fontsize, fontfamily='sans-serif')

                #Write CSV OMNI
                if data_record == True:
                    start_dtstr = self.start_dt.replace(':', '-')
                    end_dtstr = self.end_dt.replace(':', '-')
                    csv_filename = '_'.join(
                        [para, 'OMNI', start_dtstr, end_dtstr]) + '.csv'
                    # see CSA_InstrData for code
                    with open(csv_filename, 'w', newline='') as csvfile:
                        csvwriter = csv.writer(csvfile)
                        for i in np.arange(len(self.data['Epoch'].values)):
                            csvwriter.writerow([self.data['Epoch'].values[i], self.data[Bx_varname].values[i],
                                                self.data[By_varname].values[i], self.data[Bz_varname].values[i]])

                    csvfile.close()
                    print('{}, OMNI data recorded in csv'.format(para))

        #Save figure in folder sorted by start year/month/day.
        savefolder = 'OMNIPanels/{}/{}/{}/'.format(start_dtobj.year,start_dtobj.month,
                                                  start_dtobj.day)
        os.makedirs(savefolder,exist_ok=True)
        
        savename = 'OMNI'
        for para in self.paralist:
            savename = savename + '-{}'.format(para)
        savename = savename + '__{}__{}_{}.{}'.format(self.start_dt.replace(':','-'),
                                                   self.end_dt.replace(':','-'),
                                                   self.coord,saveformat)
        savepath = savefolder + savename
        
        plt.subplots_adjust(left=0.08, bottom=0.07, right=0.92, top=0.97,
                            wspace=0.0, hspace=0.08)  # 0.05 = 5% from left figure width/height
        plt.savefig(savepath, dpi=200, format=saveformat, transparent=tchoice)
        plt.close()

    def csv_record(self, para):
        #method for recording plot data as csv file
        start_dtstr = self.start_dt.strftime('%Y-%m-%dT%H-%M-%SZ')
        stop_dtstr = self.end_dt.strftime('%Y-%m-%dT%H-%M-%SZ')
        csv_filename = '_'.join(
            [para, 'OMNI', start_dtstr, stop_dtstr]) + '.csv'
        with open(csv_filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)

            for i in np.arange(len(self.data['Epoch'].values)):
                csvwriter.writerow(
                    [self.data['Epoch'].values[i], self.data[para].values[i]])

        csvfile.close()
        print('{}, OMNI data recorded in csv'.format(para))


"""
Example: Calling CDAweb manager for OMNI
"""
"""

start_dt = dt.datetime(2002,3,18,14,15,0)
start_dt = start_dt.strftime('%Y-%m-%dT%H:%M:%SZ') #Reformat to string
end_dt = dt.datetime(2002,3,18,15,15,0)
end_dt = end_dt.strftime('%Y-%m-%dT%H:%M:%SZ')
paralist = ['Bx','By','Bz']
OG = 'OMNI (Combined 1AU IP Data; Magnetic and Solar Indices)'
IT = 'Magnetic Fields (space)'

bundle = CDAWeb_Manager(start_dt,end_dt,paralist,OG,IT,coord='GSE') #Call data in GSE coordinates because no GSM data is available from CDAWeb
bundle.GSM_from_GSE('B') #Call GSM conversion manually 
bundle.OMNI_Plot(xax_view=1,label_start='A',veccomb=True,saveformat='png',
                 MP_HL=True,Walen_HL=True,data_record=True)
#xax_view: 0 = bottom-only, 1 = all invisible, 2 = all visible
"""