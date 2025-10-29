"""
Edition Date: 2025-October-27
@author: Nawapat Kaweeyanun
"""
"""
Objective: Plot Cluster's selected plasma-magnetic parameter data between selected datetimes

"""

import numpy as np
import os
import shutil
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import AutoMinorLocator
from string import ascii_uppercase as ALC
from spacepy import coordinates as coord
from spacepy.time import Ticktock
import CSA_Extract as CEx
import CSA_Organiser as Corg
import CSA_Instr_Data as CID
import csv


class CPanel(object):

    def __init__(self, start_dt, stop_dt, sclist, paralist, instrlist, saveformat,
                 new_dl=True, veccomb=False, sc_comp=False, autoplot=True,
                 MP_HL=False, Walen_HL=False, label_start='A',xax_view=0,
                 mode='light', data_record=False):
        
        self.start_dt = start_dt  #Period start datetime
        self.stop_dt = stop_dt  #Period end datetime
        self.sclist = sclist  #Input list of spacecraft
        self.paralist = paralist  #Input list of parameters
        self.instrlist = instrlist  #List of instrument to search
        self.saveformat = saveformat #Format of figure file to save 
        self.veccomb = veccomb #If true, present three vector components in one plot.
        self.MP_HL = MP_HL #If true, highlight magnetopause crossing.
        self.Walen_HL = Walen_HL #If true, highlight period used for Walen tests.
        
        #Check if instrument input is on available list. If not, remove.
        for no in np.arange(len(self.instrlist)):
            self.instrlist[no] = self.instrlist[no].upper()
            if self.instrlist[no] not in ['FGM', 'EFW', 'CIS', 'PEA', 'EDI', 'RAP']:
                print('Instrument {} not in search bank'.format(
                    self.instrlist[no]))
                self.instrlist.remove(self.instrlist[no])

        #Check that figure format is png or svg. Other formats not supported.
        if self.saveformat not in ['png', 'svg']:
            raise Exception('Figure format not supported')

        #Sort parameter by spacecraft using dictionary.
        self.paradict = {}
        for sc in self.sclist:
            self.paradict.update({sc: self.paralist})

        #Call variable sorter to obtain data IDs needed.
        self.vardict = self.VarSorter(
            self.paradict, self.instrlist, prior_sorted=False)
        print('Variable IDs Sorted')

        #Generate file download list.
        DataIDlist = self.DownloadTarget(self.vardict)
        print('Ready to Download')

        #If true, download data from Cluster Science Archive.
        #If data is already available, can set as False to save time.
        if new_dl == True:
            myurl = 'https://csa.esac.esa.int/csa-sl-tap/data'
            start_dtstr = self.start_dt.strftime(
                '%Y-%m-%dT%H:%M:%SZ')  #Remove milliseconds from datetime.
            stop_dtstr = self.stop_dt.strftime('%Y-%m-%dT%H:%M:%SZ')

            #Create intermediate data directory. Remove previous run if exists.
            if 'CSA_ToMove' in os.listdir(os.getcwd()):
                shutil.rmtree('CSA_ToMove')
            os.makedirs('CSA_ToMove', exist_ok=True)

            #For each data ID, set up parameters then request data
            for ID in DataIDlist:
                query_specs = {'RETRIEVAL_TYPE': 'product',
                               'DATASET_ID': ID,
                               'START_DATE': start_dtstr,
                               'END_DATE': stop_dtstr,
                               'DELIVERY_FORMAT': 'CDF'}
                
                try:
                    CEx.download(myurl, query_specs, 'tap_download.tgz')
                    filepath = CEx.unpack_tgz('tap_download.tgz')
                except Exception as e:
                    if 'EDI' in ID:  #For EDI can replace with lower res
                        ID = ID.split('MP')[0] + 'SPIN'
                        CEx.download(myurl, query_specs, 'tap_download.tgz')
                        filepath = CEx.unpack_tgz('tap_download.tgz')
                    else:
                        print(e)
                        print('Error downloading {}. Check if data exist'.format(ID))
                        continue
                
                #Move downloaded folder to intermediate directory
                filedir = '/'.join(filepath.split('/')[0:-1])
                #print('CSA_ToMove' + '/' + filedir)
                shutil.move(filedir, 'CSA_ToMove' + '/' + filedir)

            print('File(s) Downloaded')

            #Remove downloaded folders to reduce operating directoryclutter
            DownloadFolderList = [d for d in os.listdir(
                os.getcwd()) if 'CSA_Download' in d]
            for DLF in DownloadFolderList:
                try:
                    shutil.rmtree(DLF)
                except:  #skip undeletable files if exist
                    continue

            #Organise downloaded files into a Data Folder. Remove previous run if exists.
            if 'PanelData' in os.listdir(os.getcwd()):
                shutil.rmtree('PanelData')
            os.makedirs('PanelData', exist_ok=True)
            Corg.fileshifter('CSA_ToMove', 'PanelData')

        print('Files Organised')

        #Extract data from Data Folder. Create data dictionary
        self.DataDict = self.DictExtract('PanelData')
        print('Data extracted')
        
        #From data dictionary, rearrange data into plotting order. Create plotting dictionary.
        self.PanelDict = self.PanelArrange(
            self.paralist, self.sclist, self.DataDict)
        print('Plotting Data Obtained')

        #If set true, call plotting function
        if autoplot == True:
            self.Plot_Func(self.PanelDict, self.sclist,
                           self.MP_HL, self.Walen_HL,self.saveformat,mode,data_record)
            print('Data Plotted')

        #If set true, produce  diagram comparing spectrograms (flux plots) between
        #multiple Cluster spacecrafts (for additional analysis)
        if sc_comp == True:
            self.PanelCompDict = self.PanelCompArrange(
                self.paralist, self.sclist, self.DataDict)
            print('FluxComp Data Obtained')
            if autoplot == True:
                self.CompPlot_Func(self.PanelCompDict,
                                   self.sclist, mode, self.saveformat,data_record)
                print('FluxComp Data Plotted')

    def VarSorter(self, paradict, searchlist, prior_sorted=False):
        #Return dictionary of full variable ID names
        #Structure = Dict[para][sc]

        vardict = {}

        #Set up search net i.e., which parameters in which instruments
        #Datetimes and quality-related parameters not included
        FGM_net = ['Bvec', 'Bmag', 'SCpos']
        EFW_net = ['Evec', 'SC_Pot']
        CIS_net = ['ion_fluxU', 'ion_fluxU0', 'ion_fluxU90', 'ion_fluxU180', 'ion_fluxPA', 
                   'ion_fluxU_Sun', 'ion_fluxU_Dawn', 'ion_fluxU_Antisun', 'ion_fluxU_Dusk',
                   'Nion', 'ion_vel', 'Pion', 'Tion', 'TionBll', 'TionBll',
                   'H+_fluxU', 'H+_fluxU0', 'H+_fluxU90', 'H+_fluxU180', 'H+_fluxPA', 
                   'H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk', 
                   'NH+', 'H+_vel', 'PH+', 'TH+', 'TH+Bll', 'TH+Bp', 
                   'He+_fluxU', 'He+_fluxU0', 'He+_fluxU90', 'He+_fluxU180', 'He+_fluxPA', 
                   'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk', 
                   'NHe+', 'He+_vel', 'PHe+', 'THe+', 'THe+Bll', 'THe+_Bp', 
                   'O+_fluxU', 'O+_fluxU0', 'O+_fluxU90', 'O+_fluxU180', 'O+_fluxPA', 
                   'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk', 
                   'NO+', 'O+_vel', 'PO+', 'TO+', 'TO+Bll', 'TO+Bp']
        PEA_net = ['e_fluxU', 'e_fluxU0',
                   'e_fluxU90', 'e_fluxU180', 'e_fluxPA']
        EDI_net = ['e_vel', 'Evec']

        #If instrument not in search list, set net to empty
        if 'FGM' not in searchlist:
            FGM_net = []
        if 'EFW' not in searchlist:
            EFW_net = []
        if 'CIS' not in searchlist:
            CIS_net = []
        if 'PEA' not in searchlist:
            PEA_net = []
        if 'EDI' not in searchlist:
            EDI_net = []

        #Loop over all spacecraft inputs
        for key in paradict.keys():
            vardict[key] = []

            #Loop over parameter list
            for para in paradict[key]:

                para_priority = {}

                #Some parameters can be obtained from more than one instruments. If desire, set priority
                #False: EFW first then EDI
                if para in ['Evec'] and prior_sorted == False:

                    C = input(
                        'Priority for {}? Unless type "no", default to yes'.format(para))

                    if C.lower() == 'no':
                        para_priority.update({para: False})
                        self.Eprior = False
                    else:
                        para_priority.update({para: True})
                        self.Eprior = True

                #If already sorted once, no need to do again
                elif para in ['Evec'] and prior_sorted == True:
                    para_priority.update({para: self.Eprior})

                #Obtain data ID for the parameter.
                #Use spin resolution data over full where relevant.
                #Use high resolution flux data where relevant
                if para in FGM_net:
                    if para == 'Bvec':
                        paraID = 'B_vec_xyz_gse__{}_CP_FGM_SPIN'.format(
                            key.upper())
                    elif para == 'Bmag':
                        paraID = 'B_mag__{}_CP_FGM_SPIN'.format(key.upper())
                    elif para == 'SCpos':
                        paraID = 'sc_pos_xyz_gse__{}_CP_FGM_SPIN'.format(
                            key.upper())
                elif para in EFW_net:
                    if para == 'Evec' and para_priority['Evec'] == True:
                        paraID = 'E_Vec_xyz_GSE__{}_CP_EFW_L3_E3D_GSE'.format(
                            key.upper())
                    elif para == 'Evec' and para_priority['Evec'] == False:
                        paraID = 'E_xyz_gse__{}_CP_EDI_MP'.format(
                            key.upper())
                    elif para == 'SC_Pot':
                        paraID = 'Spacecraft_potential__{}_CP_EFW_L3_P'.format(
                            key.upper())
                elif para in CIS_net:
                    if 'ion' in para:
                        if para in ['ion_fluxU', 'ion_fluxU0', 'ion_fluxU90', 'ion_fluxU180', 'ion_fluxPA']:
                            paraID = 'Differential_Particle_fluxU__{}_CP_CIS-HIA_PAD_HS_MAG_IONS_PF'.format(
                                key.upper())
                        elif para in ['ion_fluxU_Sun', 'ion_fluxU_Dawn', 'ion_fluxU_Antisun', 'ion_fluxU_Dusk']:
                            paraID = '3d_ions__{}_CP_CIS-HIA_HS_MAG_IONS_PF'.format(
                                key.upper())
                        elif para == 'Nion':
                            paraID = 'density__{}_CP_CIS-HIA_ONBOARD_MOMENTS'.format(
                                key.upper())
                        elif para == 'ion_vel':
                            paraID = 'velocity__{}_CP_CIS-HIA_ONBOARD_MOMENTS'.format(
                                key.upper())
                        elif para == 'Pion':
                            paraID = 'pressure__{}_CP_CIS-HIA_ONBOARD_MOMENTS'.format(
                                key.upper())
                        elif para == 'Tion':
                            paraID = 'temperature__{}_CP_CIS-HIA_ONBOARD_MOMENTS'.format(
                                key.upper())
                        elif para == 'TionBll':
                            paraID = 'temp_par__{}_CP_CIS-HIA_ONBOARD_MOMENTS'.format(
                                key.upper())
                        elif para == 'TionBll':
                            paraID = 'temp_perp__{}_CP_CIS-HIA_ONBOARD_MOMENTS'.format(
                                key.upper())

                    elif 'H+' in para:
                        if para in ['H+_fluxU', 'H+_fluxU0', 'H+_fluxU90', 'H+_fluxU180', 'H+_fluxPA']:
                            paraID = 'Differential_Particle_fluxU__{}_CP_CIS-CODIF_PAD_HS_H1_PF'.format(
                                key.upper())
                        elif para in ['H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk']:
                            paraID = '3d_ions__{}_CP_CIS-CODIF_HS_H1_PF'.format(
                                key.upper())
                        elif para == 'NH+':
                            paraID = 'density__{}_CP_CIS-CODIF_HS_H1_MOMENTS'.format(
                                key.upper())
                        elif para == 'H+_vel':
                            paraID = 'velocity__{}_CP_CIS-CODIF_HS_H1_MOMENTS'.format(
                                key.upper())
                        elif para == 'PH+':
                            paraID = 'pressure__{}_CP_CIS-CODIF_HS_H1_MOMENTS'.format(
                                key.upper())
                        elif para == 'TH+':
                            paraID = 'T__{}_CP_CIS-CODIF_HS_H1_MOMENTS'.format(
                                key.upper())
                        elif para == 'TH+Bll':
                            paraID == 'T_par__{}_CP_CIS-CODIF_HS_H1_MOMENTS'.format(
                                key.upper())
                        elif para == 'TH+Bp':
                            paraID == 'T_perp__{}_CP_CIS-CODIF_HS_H1_MOMENTS'.format(
                                key.upper())

                    elif 'He+' in para:
                        if para in ['He+_fluxU', 'He+_fluxU0', 'He+_fluxU90', 'He+_fluxU180', 'He+_fluxPA']:
                            paraID = 'Differential_Particle_fluxU__{}_CP_CIS-CODIF_PAD_HS_He1_PF'.format(
                                key.upper())        
                        elif para in ['He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk']:
                            paraID = '3d_ions__{}_CP_CIS-CODIF_HS_HE1_PF'.format(
                                key.upper())
                        elif para == 'NHe+':
                            paraID = 'density__{}_CP_CIS-CODIF_HS_He1_MOMENTS'.format(
                                key.upper())
                        elif para == 'He+_vel':
                            paraID = 'velocity__{}_CP_CIS-CODIF_HS_He1_MOMENTS'.format(
                                key.upper())
                        elif para == 'PHe+':
                            paraID = 'pressure__{}_CP_CIS-CODIF_HS_He1_MOMENTS'.format(
                                key.upper())
                        elif para == 'THe+':
                            paraID = 'T__{}_CP_CIS-CODIF_HS_He1_MOMENTS'.format(
                                key.upper())
                        elif para == 'THe+Bll':
                            paraID == 'T_par__{}_CP_CIS-CODIF_HS_H1_MOMENTS'.format(
                                key.upper())
                        elif para == 'THe+Bp':
                            paraID == 'T_perp__{}_CP_CIS-CODIF_HS_H1_MOMENTS'.format(
                                key.upper())

                    elif 'O+' in para:
                        if para in ['O+_fluxU', 'O+_fluxU0', 'O+_fluxU90', 'O+_fluxU180', 'O+_fluxPA']:
                            paraID = 'Differential_Particle_fluxU__{}_CP_CIS-CODIF_PAD_HS_O1_PF'.format(
                                key.upper())
                        elif para in ['O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk']:
                            paraID = '3d_ions__{}_CP_CIS-CODIF_HS_O1_PF'.format(
                                key.upper())
                        elif para == 'NO+':
                            paraID = 'density__{}_CP_CIS-CODIF_HS_O1_MOMENTS'.format(
                                key.upper())
                        elif para == 'O+_vel':
                            paraID = 'velocity__{}_CP_CIS-CODIF_HS_O1_MOMENTS'.format(
                                key.upper())
                        elif para == 'PO+':
                            paraID = 'pressure__{}_CP_CIS-CODIF_HS_O1_MOMENTS'.format(
                                key.upper())
                        elif para == 'TO+':
                            paraID = 'T__{}_CP_CIS-CODIF_HS_O1_MOMENTS'.format(
                                key.upper())
                        elif para == 'TO+Bll':
                            paraID == 'T_par__{}_CP_CIS-CODIF_HS_H1_MOMENTS'.format(
                                key.upper())
                        elif para == 'TO+Bp':
                            paraID == 'T_perp__{}_CP_CIS-CODIF_HS_H1_MOMENTS'.format(
                                key.upper())

                elif para in PEA_net:
                    if para in ['e_fluxU', 'e_fluxU0', 'e_fluxU90', 'e_fluxU180', 'e_fluxPA']:
                        paraID = 'Data__{}_CP_PEA_PITCH_SPIN_DEFlux'.format(
                            key.upper())

                elif para in EDI_net:
                    #Note: derived, not calibrated data
                    if para == 'e_vel':
                        paraID = 'V_ed_xyz_gse__{}_CP_EDI_MP'.format(
                            key.upper())
                    elif para == 'Evec':
                        paraID = 'E_xyz_gse__{}_CP_EDI_MP'.format(key.upper())

                else:
                    print('Parameter {} not found in given search bank'.format(para))
                    paraID = 'N/A'

                vardict[key].append(paraID)

        return vardict

    def DownloadTarget(self, vardict):
        #Generate a list of files to download from given key variable dictionary
        downloadlist = []
        for key in vardict.keys():
            for para in vardict[key]:
                FileID = para.split('__')[1]
                downloadlist.append(FileID)

        #Remove identical elements. Order is scrambled but doesn't matter.
        downloadlist = list(set(downloadlist))

        return downloadlist

    def DictExtract(self, sourcefolder):
        #Extract data from all relevant files
        #Final dictionary has format Dict[DataID][datetime-period]

        Datadict = {}

        #Loop over every possible file in Data Folder
        IDlist = os.listdir(sourcefolder)
        IDlist.sort()

        for ID in IDlist:

            IDdict = {}

            for yr in os.listdir(sourcefolder + '/' + ID):
                for mo in os.listdir(sourcefolder + '/{}/{}'.format(ID, yr)):
                    filelist = os.listdir(
                        sourcefolder + '/{}/{}/{}'.format(ID, yr, mo))

                    for file in filelist:
                        filepath = sourcefolder + \
                            '/{}/{}/{}/{}'.format(ID, yr, mo, file)
                        DataSet = CID.Cdata(filepath)

                        filedict = {}

                        #Loop over all data variables and append non-functions
                        for att in dir(DataSet):

                            #Set up the string for determining callability
                            str_to_exec = """if callable(DataSet.%s) == False and '__' not in '%s':
                                value = DataSet.%s
                                filedict.update({'%s':value})
                            """
                            #Execute the string
                            exec(str_to_exec % (att, att, att, att))

                        period_str = filedict['periodname']
                        IDdict.update({period_str: filedict})

            #Update data dictionary
            Datadict.update({ID: IDdict})

        return Datadict

    def PanelArrange(self, paralist, sclist, DataDict):
        #Set up by panel plotting dictionarydictionary.
        #Format Dict[para][sc] if data can be plotted for multiple spacecrafts in a panel (Multiple)
        #Format Dict[para] if data can be plotted for single spacecraft only (Single, flux data)

        PanelDict = {}

        Multiple = ['Bvec', 'Bmag', 'SCpos', 'Evec', 'Nion', 'ion_vel', 'Pion',
                    'Tion', 'TionBll', 'TionBll', 'NH+', 'H+_vel',
                    'PH+', 'TH+', 'TH+Bll', 'TH+Bp', 'NHe+',
                    'He+_vel', 'PHe+', 'THe+', 'THe+Bll', 'THe+Bp',
                    'NO+', 'O+_vel', 'PO+', 'TO+', 'TO+Bll', 'TO+Bp',
                    'e_vel', 'Evec', 'SC_Pot']
        Single = ['ion_fluxU', 'ion_fluxU0', 'ion_fluxU90', 'ion_fluxU180', 'ion_fluxPA', 
                  'ion_fluxU_Sun', 'ion_fluxU_Dawn', 'ion_fluxU_Antisun', 'ion_fluxU_Dusk', 
                  'H+_fluxU', 'H+_fluxU0', 'H+_fluxU90', 'H+_fluxU180', 'H+_fluxPA',
                  'H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk', 
                  'He+_fluxU', 'He+_fluxU0', 'He+_fluxU90', 'He+_fluxU180', 'He+_fluxPA', 
                  'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk', 
                  'O+_fluxU', 'O+_fluxU0', 'O+_fluxU90', 'O+_fluxU180', 'O+_fluxPA', 
                  'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk', 
                  'e_fluxU', 'e_fluxU0', 'e_fluxU90', 'e_fluxU180', 'e_fluxPA']
        
        #Special list for parameters requiring 3D particle distributions
        Para3D = ['ion_fluxU_Sun', 'ion_fluxU_Dawn', 'ion_fluxU_Antisun', 'ion_fluxU_Dusk', 
                    'H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk', 
                    'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk', 
                    'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk']

        #Set up final dictionary
        for para in paralist:
            if para not in Multiple and para not in Single: #remove parameters not in either category
                print('Parameter {} not covered. Removed.'.format(para))
                paralist.remove(para)  #only work inside this function
                continue
            if para == 'Bvec' and self.veccomb == False:  #separate vector components if chosen so
                PanelDict.update({'Bx': {}})
                PanelDict.update({'By': {}})
                PanelDict.update({'Bz': {}})
                for sc in sclist:
                    PanelDict['Bx'].update({sc: []})
                    PanelDict['By'].update({sc: []})
                    PanelDict['Bz'].update({sc: []})
            elif para == 'Evec' and self.veccomb == False:
                PanelDict.update({'E': {}})
                PanelDict.update({'Ex': {}})
                PanelDict.update({'Ey': {}})
                PanelDict.update({'Ez': {}})
                for sc in sclist:
                    PanelDict['E'].update({sc: []})
                    PanelDict['Ex'].update({sc: []})
                    PanelDict['Ey'].update({sc: []})
                    PanelDict['Ez'].update({sc: []})
            elif para == 'SCpos' and self.veccomb == False:
                PanelDict.update({'Rgse': {}})
                PanelDict.update({'xgse': {}})
                PanelDict.update({'ygse': {}})
                PanelDict.update({'zgse': {}})
                for sc in sclist:
                    PanelDict['Rgse'].update({sc: []})
                    PanelDict['xgse'].update({sc: []})
                    PanelDict['ygse'].update({sc: []})
                    PanelDict['zgse'].update({sc: []})
            elif para == 'ion_vel' and self.veccomb == False:
                PanelDict.update({'vi': {}})
                PanelDict.update({'vxi': {}})
                PanelDict.update({'vyi': {}})
                PanelDict.update({'vzi': {}})
                for sc in sclist:
                    PanelDict['vi'].update({sc: []})
                    PanelDict['vxi'].update({sc: []})
                    PanelDict['vyi'].update({sc: []})
                    PanelDict['vzi'].update({sc: []})
            elif para == 'H+_vel' and self.veccomb == False:
                PanelDict.update({'vH': {}})
                PanelDict.update({'vxH': {}})
                PanelDict.update({'vyH': {}})
                PanelDict.update({'vzH': {}})
                for sc in sclist:
                    PanelDict['vH'].update({sc: []})
                    PanelDict['vxH'].update({sc: []})
                    PanelDict['vyH'].update({sc: []})
                    PanelDict['vzH'].update({sc: []})
            elif para == 'He+_vel' and self.veccomb == False:
                PanelDict.update({'vHe': {}})
                PanelDict.update({'vxHe': {}})
                PanelDict.update({'vyHe': {}})
                PanelDict.update({'vzHe': {}})
                for sc in sclist:
                    PanelDict['vHe'].update({sc: []})
                    PanelDict['vxHe'].update({sc: []})
                    PanelDict['vyHe'].update({sc: []})
                    PanelDict['vzHe'].update({sc: []})
            elif para == 'O+_vel' and self.veccomb == False:
                PanelDict.update({'vO': {}})
                PanelDict.update({'vxO': {}})
                PanelDict.update({'vyO': {}})
                PanelDict.update({'vzO': {}})
                for sc in sclist:
                    PanelDict['vO'].update({sc: []})
                    PanelDict['vxO'].update({sc: []})
                    PanelDict['vyO'].update({sc: []})
                    PanelDict['vzO'].update({sc: []})
            elif para == 'e_vel' and self.veccomb == False:
                PanelDict.update({'ve': {}})
                PanelDict.update({'vxe': {}})
                PanelDict.update({'vye': {}})
                PanelDict.update({'vze': {}})
                for sc in sclist:
                    PanelDict['ve'].update({sc: []})
                    PanelDict['vxe'].update({sc: []})
                    PanelDict['vye'].update({sc: []})
                    PanelDict['vze'].update({sc: []})
            else:
                PanelDict.update({para: {}})
                for sc in sclist:
                    PanelDict[para].update({sc: []})

        #For each parameter, search data dictionary for correct ID, then create
        #'combined' array containing all data to be plotted.
        for para in paralist:

            if para in Multiple: #For data other than particle flux

                for key in DataDict.keys():

                    for sc in sclist:
                        convert_in = {sc: [para]}
                        convert_out = self.VarSorter(
                            convert_in, self.instrlist, prior_sorted=True)

                        ID = convert_out[sc][0].split(
                            '__')[1]

                        if ID == key:

                            print('Key {} ID Matched'.format(ID))

                            #set up array for datetime and data
                            dt_arr = np.empty([1])
                            if para in ['Bvec', 'Evec', 'SCpos', 'ion_vel', 'H+_vel', 'He+_vel,', 'O+_vel', 'e_vel']:
                                data_arr = np.empty([1, 3]) #vector
                            else: 
                                data_arr = np.empty([1])  #scalar
                            
                            #Loop over all downloaded file(s)
                            for period in DataDict[ID].keys():
                                dtime = DataDict[ID][period]['datetime']

                                if para == 'Bvec':
                                    pdata = DataDict[ID][period]['Bvec']
                                elif para == 'Bmag':
                                    pdata = DataDict[ID][period]['Bmag']
                                elif para == 'SCpos':
                                    pdata = DataDict[ID][period]['SC_gse']
                                    pdata = pdata/6371.2  # convert to RE
                                elif para == 'Evec':
                                    pdata = DataDict[ID][period]['Evec']
                                elif para in ['Nion', 'NH+', 'NHe+', 'NO+']:
                                    pdata = DataDict[ID][period]['density']
                                elif para == 'ion_vel':
                                    pdata = DataDict[ID][period]['vel_gse']
                                elif para in ['H+_vel', 'He+_vel', 'O+_vel', 'e_vel']:
                                    pdata = DataDict[ID][period]['vel']
                                elif para in ['Pion', 'PH+', 'PHe+', 'PO+']:
                                    pdata = DataDict[ID][period]['pres']
                                elif para in ['Tion', 'TH+', 'THe+', 'TO+']:
                                    pdata = DataDict[ID][period]['temp']
                                    pdata = pdata*10**6*8.625*10**-8  # convert from MK to keV
                                elif para in ['TionBll', 'TH+Bll', 'THe+Bll', 'TO+Bll']:
                                    pdata = DataDict[ID][period]['temp_Bll']
                                elif para in ['TionBll', 'TH+Bp', 'THe+Bp', 'TO+Bp']:
                                    pdata = DataDict[ID][period]['temp_Bp']
                                elif para == 'SC_Pot':
                                    pdata = DataDict[ID][period]['SC_Pot']
                                    #For datetimes with fill values, use average of nearest 10 data points
                                    #10 because can have long stretch with fill values
                                    fi_list = np.where(pdata == -1*10**9)[0]
                                    for fi in fi_list:
                                        around = pdata[fi-5:fi+5]
                                        pdata[fi] = np.mean(
                                            around, where=around != -1*10**9)

                                #Append arrays
                                dt_arr = np.concatenate((dt_arr, dtime))

                                if len(np.shape(data_arr)) == 1:
                                    #Where there is no data in 1D array (or just fill value), scalar para produces empty array[]
                                    if np.size(pdata) == 0 or np.size(pdata) == 1:
                                        data_arr = np.concatenate(
                                            (data_arr, np.array([np.nan])), axis=0)
                                    else:
                                        data_arr = np.concatenate(
                                            (data_arr, pdata), axis=0)
                                elif len(np.shape(data_arr)) == 2:
                                    #Where there is no data in 2D array, vector para produnces filler 1D array of three -1*10E31 elements
                                    if len(np.shape(pdata)) == 1:
                                        nan_arr = np.zeros((1, 3))
                                        nan_arr[:] = np.nan
                                        data_arr = np.concatenate(
                                            (data_arr, nan_arr), axis=0)
                                    else:
                                        data_arr = np.concatenate(
                                            (data_arr, pdata), axis=0)
                            
                            #Delete randomised initialising elements
                            dt_arr = np.delete(dt_arr, (0), axis=0)
                            data_arr = np.delete(data_arr, (0), axis=0)
                            
                            
                            #GSM conversion (works for both 1D and 2D data as long as dtime is along row)
                            if para in ['Bvec','Evec','ion_vel','H+_vel','He+_vel','O+_vel','e_vel','SCpos']:
                                data_arr = self.GSM_from_GSE(dt_arr, data_arr)
                            

                            #Update data into dictionary. If clause accounts for vector data.
                            if len(np.shape(data_arr)) == 2 and self.veccomb == False:
                                data_arr_x = data_arr[:, 0]  # first column
                                data_arr_y = data_arr[:, 1]
                                data_arr_z = data_arr[:, 2]

                                comb_arr_x = np.stack(
                                    (dt_arr, data_arr_x)).transpose()
                                comb_arr_y = np.stack(
                                    (dt_arr, data_arr_y)).transpose()
                                comb_arr_z = np.stack(
                                    (dt_arr, data_arr_z)).transpose()

                                if para != 'Bvec':
                                    data_arr = np.sqrt(
                                        data_arr_x**2 + data_arr_y**2 + data_arr_z**2)
                                    comb_arr = np.stack(
                                        (dt_arr, data_arr)).transpose()

                                if para == 'Bvec':
                                    PanelDict['Bx'][sc] = comb_arr_x
                                    PanelDict['By'][sc] = comb_arr_y
                                    PanelDict['Bz'][sc] = comb_arr_z
                                elif para == 'Evec':
                                    PanelDict['E'][sc] = comb_arr
                                    PanelDict['Ex'][sc] = comb_arr_x
                                    PanelDict['Ey'][sc] = comb_arr_y
                                    PanelDict['Ez'][sc] = comb_arr_z
                                elif para == 'SCpos':
                                    PanelDict['Rgse'][sc] = comb_arr
                                    PanelDict['xgse'][sc] = comb_arr_x
                                    PanelDict['ygse'][sc] = comb_arr_y
                                    PanelDict['zgse'][sc] = comb_arr_z
                                elif para == 'ion_vel':
                                    PanelDict['vi'][sc] = comb_arr
                                    PanelDict['vxi'][sc] = comb_arr_x
                                    PanelDict['vyi'][sc] = comb_arr_y
                                    PanelDict['vzi'][sc] = comb_arr_z
                                elif para == 'H+_vel':
                                    PanelDict['vH'][sc] = comb_arr
                                    PanelDict['vxH'][sc] = comb_arr_x
                                    PanelDict['vyH'][sc] = comb_arr_y
                                    PanelDict['vzH'][sc] = comb_arr_z
                                elif para == 'He+_vel':
                                    PanelDict['vHe'][sc] = comb_arr
                                    PanelDict['vxHe'][sc] = comb_arr_x
                                    PanelDict['vyHe'][sc] = comb_arr_y
                                    PanelDict['vzHe'][sc] = comb_arr_z
                                elif para == 'O+_vel':
                                    PanelDict['vO'][sc] = comb_arr
                                    PanelDict['vxO'][sc] = comb_arr_x
                                    PanelDict['vyO'][sc] = comb_arr_y
                                    PanelDict['vzO'][sc] = comb_arr_z
                                elif para == 'e_vel':
                                    PanelDict['ve'][sc] = comb_arr
                                    PanelDict['vxe'][sc] = comb_arr_x
                                    PanelDict['vye'][sc] = comb_arr_y
                                    PanelDict['vze'][sc] = comb_arr_z

                            else:
                                if len(np.shape(data_arr)) == 1: #if 1D data array
                                    comb_arr = np.stack(
                                        (dt_arr, data_arr)).transpose()
                                    PanelDict[para][sc] = comb_arr
                                elif len(np.shape(data_arr)) == 2: #if 2D data array
                                    dt_arr = dt_arr[:, None]
                                    comb_arr = np.concatenate(
                                        (dt_arr, data_arr), axis=1)
                                    PanelDict[para][sc] = comb_arr
                        else:
                            continue

            elif para in Single: # for particle flux data
                
                #When multiple spacecrafts are entered, only one spacecraft data can be published
                #in a panel. The next section utilises a while loop to iterate until the first valid
                #i.e., has data, spacecraft is found. Once found, the switch is set to false and the
                #loop terminates.
                
                #The code constructs one combined 2D flux data array containing all plottable information
                #for each parameter
            
                for key in DataDict.keys():
                    no = 0
                    switch = True

                    while switch == True and no < len(sclist):
                        sc = sclist[no]
                        convert_in = {sc: [para]}
                        convert_out = self.VarSorter(
                            convert_in, self.instrlist, prior_sorted=True)

                        ID = convert_out[sc][0].split(
                            '__')[1]

                        if ID != key:
                            no = no+1
                            continue
                        else:
                            print('Key {} ID Matched'.format(ID))

                            dt_arr = np.empty([1])
                            energy_arr = np.empty([1])
                            pa_arr = np.empty([1])

                            if 'CIS' in ID:
                                if 'U' in para:
                                    # 16 = pa, 31 = energy
                                    data_arr = np.empty([1, 31])
                                elif 'PA' in para:
                                    data_arr = np.empty([1, 16])
                            elif 'PEA' in ID:
                                if 'U' in para:
                                    #12 = pa, 44 = energy
                                    data_arr = np.empty([1, 44])
                                elif 'PA' in para:
                                    data_arr = np.empty([1, 12])
                            
                            #Append all available data together from file(s)
                            for period in DataDict[ID].keys():
                                
                                dtime = DataDict[ID][period]['datetime']
                                #Stop uploading corrupt datetime
                                if isinstance(dtime[0], dt.datetime) == False:
                                    print(
                                        'Append Stop Due to Datetime Issue for ID {} and Period {}'.format(ID, period))
                                    continue
                                
                                #Remove fill values for energy array and convert to keV except in case of pitch angle plots
                                energy = DataDict[ID][period]['energy']
                                if 'PEA' in ID:
                                    if len(np.shape(energy)) != 2:  #Filled values result in 1D arrays
                                        continue
                                    energy = energy[0, :]/1000 
                                else:
                                    if 'fluxU' in para:
                                        energy = energy/1000
                                    else:
                                        energy = energy
                                
                                #Set up angle arrays in degrees
                                if para in Para3D: #3D distributions
                                    theta = DataDict[ID][period]['theta']
                                    phi = DataDict[ID][period]['phi']
                                else: #Pitch Angle distributions
                                    pa = DataDict[ID][period]['pa']
                                
                                #Make a copy of differential energy flux so global value are not changed if array is manipulated more than once
                                #(Occur in Fig S2)
                                pdata_init = np.copy(DataDict[ID][period]['diffflux'])
                                
                                #Remove filled values that change array dimensions
                                if para in Para3D:
                                    if len(np.shape(pdata_init)) != 4: # remove 3D filled values
                                        continue
                                else:
                                    if len(np.shape(pdata_init)) != 3: # remove 2D filled values
                                        continue

                                #For CIS, convert particle flux to energy flux
                                if 'CIS' in ID:
                                    if para in Para3D:
                                        for i in np.arange(len(dtime)):
                                            for j in np.arange(len(theta)):
                                                for k in np.arange(len(phi)):
                                                    pdata_init[i,j,k,:] = pdata_init[i,j,k,:]*energy #In keV/cm2 s sr keV
                                    
                                    else:
                                        for i in np.arange(len(dtime)):
                                            for j in np.arange(len(pa)):
                                                if 'fluxU' in para:
                                                    pdata_init[i, j,
                                                               :] = pdata_init[i, j, :]*energy
                                                else:
                                                    #Must convert to keV when conversion was not made above.
                                                    pdata_init[i, j, :] = pdata_init[i,
                                                                                     j, :]*energy/1000

                                if 'fluxU' in para and para not in Para3D:  #Flux by energy from pitch angle distributions
                                    pdata = np.zeros((len(dtime), len(energy)))
                                    
                                    #Average over all pitch angles for general flux data.
                                    #If pitch angle is specified, pick relevant channels.
                                    #PA = 90 data is averaged over two middle channels

                                    for i in np.arange(len(dtime)):
                                        for j in np.arange(len(energy)):

                                            if 'U0' in para:
                                                pdata[i, j] = pdata_init[i, 0, j]
                                            elif 'U90' in para:
                                                if 'CIS' in ID:
                                                    pdata_pa = pdata_init[i, 7:9, j]
                                                elif 'PEA' in ID:
                                                    pdata_pa = pdata_init[i, 5:7, j]
                                                pdata[i, j] = np.sum(
                                                    pdata_pa)/len(pdata_pa)
                                            elif 'U180' in para:
                                                pdata[i, j] = pdata_init[i, -1, j]
                                            else:

                                                pdata_pa = pdata_init[i, :, j]
                                                pdata[i, j] = np.sum(
                                                    pdata_pa)/len(pdata_pa)

                                    dt_arr = np.concatenate((dt_arr, dtime))
                                    energy_arr = np.concatenate(
                                        (energy_arr, energy))
                                    data_arr = np.concatenate(
                                        (data_arr, pdata), axis=0)

                                elif 'fluxPA' in para:  #Flux by pitch angle from pitch angle distributions
                                    pdata = np.zeros((len(dtime), len(pa)))

                                    for i in np.arange(len(dtime)):
                                        for j in np.arange(len(pa)):
                                            pdata_U = pdata_init[i, j, :]

                                            if 'SC_Pot' in paralist and 'e_flux' in para:  
                                                #Remove photo-electrons by using spacecraft potential data from EFW if spacecraft potential is entered as input.
                                                #Spacecraft potential is dynamic and varies with datetime.
                                                
                                                PotID = '{}_CP_EFW_L3_P'.format(
                                                    sc.upper())
                                                Pot_Period = None
                                                
                                                #Match time period to obtain correct spacecraft potential data
                                                for c in DataDict[PotID].keys():
                                                    PP_startstr = c.split(
                                                        '___')[0]
                                                    PP_start = dt.datetime.strptime(
                                                        PP_startstr, '%Y-%m-%dT%H-%M-%S')
                                                    PP_endstr = c.split('___')[
                                                        1]
                                                    PP_end = dt.datetime.strptime(
                                                        PP_endstr, '%Y-%m-%dT%H-%M-%S')

                                                    if dtime[i] > PP_start and dtime[i] < PP_end:
                                                        Pot_Period = c
                                                if Pot_Period == None:
                                                    print(
                                                        'Photo-Electron Energy Threshold not Found')
                                                    pdata[i, j] = np.sum(
                                                        pdata_U)/len(pdata_U)
                                                    continue

                                                #Find photo-electron energy at closest datetime (good enough except very near magnetopause)
                                                Pot_dtime = DataDict[PotID][Pot_Period]['datetime']
                                                diff_arr = abs(
                                                    dtime[i]-Pot_dtime)
                                                match_ind = np.where(
                                                    diff_arr == min(diff_arr))
                                                match_Pot = DataDict[PotID][Pot_Period]['SC_Pot'][match_ind]

                                                photoe = -1*match_Pot/1000  # In keV
                                                photoe = photoe*1.30  #Add 30% threshold for boom-body potential difference
                                                #If more than one match, pick max E
                                                if len(photoe) > 1:
                                                    photoe = max(photoe)
                                                pdata_Unew = pdata_U[energy > photoe]

                                                if len(pdata_Unew) != 0:
                                                    pdata[i, j] = np.sum(
                                                        pdata_Unew)/len(pdata_U)
                                                elif len(pdata_Unew) == 0:
                                                    pdata[i, j] = np.sum(
                                                        pdata_U)/len(pdata_U)
                                            else:
                                                pdata[i, j] = np.sum(
                                                    pdata_U)/len(pdata_U)

                                    dt_arr = np.concatenate((dt_arr, dtime))
                                    pa_arr = np.concatenate((pa_arr, pa))
                                    data_arr = np.concatenate(
                                        (data_arr, pdata), axis=0)  # should get dt x pa
                                
                                elif para in Para3D: #For by energy plots requiring 3D particle distributions
                                    pdata = np.zeros((len(dtime), len(energy)))
                                    for i in np.arange(len(dtime)):
                                        for j in np.arange(len(energy)):
                                            
                                            #Pick azimuthal angle (phi) to sum according to ISR2 coordinates
                                            #Note: angle data is in direction of arrival so dawn and dusk directions are reflected from coordinate definitions 
                                            #Obtain 2D in [theta, phi] dimension
                                            if 'fluxU_Sun' in para:
                                                pdata_th = pdata_init[i,:,2:10,j] #5:7 direct S
                                            elif 'fluxU_Dawn' in para:
                                                pdata_th = pdata_init[i,:,6:14,j] #9:11 direct dawn
                                            elif 'fluxU_Antisun' in para:
                                                pdata_thA = pdata_init[i,:,10:,j] #13:15 direct AS
                                                pdata_thB = pdata_init[i,:,0:2,j] #there is a disconnect between relevant channels
                                                pdata_th = np.concatenate((pdata_thA,pdata_thB),axis=1)
                                            elif 'fluxU_Dusk' in para:
                                                pdata_thA = pdata_init[i,:,-2:,j] #2:4 direct dusk
                                                pdata_thB = pdata_init[i,:,0:6,j]
                                                pdata_th = np.concatenate((pdata_thA,pdata_thB),axis=1)
                                            
                                            #Average over grouped azimuthal angles.
                                            pdata_th = np.sum(pdata_th,axis=1)/len(pdata_th[0,:]) #average along horizontal phi axis
                                            
                                            #Then average over elevation angles
                                            pdata[i,j] = np.sum(pdata_th)/len(pdata_th)
                                            
                                    dt_arr = np.concatenate((dt_arr, dtime))
                                    energy_arr = np.concatenate(
                                        (energy_arr, energy))
                                    data_arr = np.concatenate(
                                        (data_arr, pdata), axis=0)
                                
                            #Remove initialising first element
                            dt_arr = np.delete(dt_arr, (0), axis=0)
                            energy_arr = np.delete(energy_arr, (0), axis=0)
                            pa_arr = np.delete(pa_arr, (0), axis=0)
                            data_arr = np.delete(data_arr, (0), axis=0)

                            #Fill the end of either energy or datetime arrays with NaN to be stackable
                            #Fill the array with shorter length. If datetime is filled, also must extend data array
                            if 'fluxU' in para:

                                #Potential addition for when energy array is longer than dt array
                                if len(dt_arr) >= len(energy_arr):
                                    nan_extra = np.nan * \
                                        np.zeros((len(dt_arr)-len(energy_arr)))
                                    energy_arr = np.concatenate(
                                        (energy_arr, nan_extra))
                                else:
                                    nan_extra = np.nan * \
                                        np.zeros((len(energy_arr)-len(dt_arr)))
                                    nan_extra2D = np.nan * \
                                        np.zeros(
                                            (len(energy_arr)-len(dt_arr), len(data_arr[0, :])))

                                    dt_arr = np.concatenate(
                                        (dt_arr, nan_extra))
                                    data_arr = np.concatenate(
                                        (data_arr, nan_extra2D), axis=0)

                                    #Need to add nan_rows to data_arr as well because data_arr always have same length as dt_arr and now energy_arr is longest
                                    #Cannot cut off energy_arr

                                comb1 = np.stack(
                                    (dt_arr, energy_arr)).transpose() #Make 2D
                                comb2 = np.concatenate(
                                    (comb1, data_arr), axis=1) #Stitch along column

                            elif 'fluxPA' in para:

                                if len(dt_arr) >= len(pa_arr):
                                    nan_extra = np.nan * \
                                        np.zeros((len(dt_arr)-len(pa_arr)))
                                    pa_arr = np.concatenate(
                                        (pa_arr, nan_extra))
                                else:
                                    nan_extra = np.nan * \
                                        np.zeros((len(pa_arr)-len(dt_arr)))
                                    nan_extra2D = np.nan * \
                                        np.zeros(
                                            (len(pa_arr)-len(dt_arr), len(data_arr[0, :])))

                                    dt_arr = np.concatenate(
                                        (dt_arr, nan_extra))
                                    data_arr = np.concatenate(
                                        (data_arr, nan_extra2D), axis=0)

                                comb1 = np.stack(
                                    (dt_arr, pa_arr)).transpose()
                                comb2 = np.concatenate(
                                    (comb1, data_arr), axis=1)

                            if np.size(data_arr) == 0: #If there is no element in arr, shift to next spacecraft
                                no = no+1
                            else: #If there is data, append to dictionary and stop the while loop
                                PanelDict[para][sc] = comb2
                                switch = False

        return PanelDict

    def Plot_Func(self, PanelDict, sclist, MP_HL, Walen_HL,saveformat,xax_view=0,mode='light',data_record=False):
        #Plot data from Plotting Dictionary

        #Define multiple & single spacecraft parameters (including vector components
        Multiple = ['Bmag', 'Bx', 'By', 'Bz', 'SCpos', 'Rgse', 'xgse', 'ygse', 'zgse', 'E', 'Ex', 'Ey', 'Ez', 'vi', 'vxi', 'vyi', 'vzi',
                    'vH', 'vxH', 'vyH', 'vzH', 'vHe', 'vxHe', 'vyHe', 'vzHe', 'vO', 'vxO', 'vyO', 'vzO', 've', 'vxe', 'vye', 'vze',
                    'Nion', 'Pion', 'Tion', 'TionBll', 'TionBll', 'NH+', 'PH+', 'TH+', 'TH+Bll', 'TH+Bp', 'NHe+', 'PHe+', 'THe+',
                    'THe+Bll', 'THe+Bp', 'NO+', 'PO+', 'TO+', 'TO+Bll', 'TO+Bp', 'Bvec', 'Evec', 'ion_vel', 'H+_vel',
                    'He+_vel', 'O+_vel', 'e_vel']
    
        Single = ['ion_fluxU', 'ion_fluxU0', 'ion_fluxU90', 'ion_fluxU180', 'ion_fluxPA', 'ion_fluxU_Sun', 'ion_fluxU_Dawn', 
                  'ion_fluxU_Antisun', 'ion_fluxU_Dusk', 'H+_fluxU', 'H+_fluxU0', 'H+_fluxU90', 'H+_fluxU180', 'H+_fluxPA', 
                  'H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk', 'He+_fluxU', 'He+_fluxU0', 'He+_fluxU90',
                  'He+_fluxU180', 'He+_fluxPA', 'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk', 'O+_fluxU', 
                  'O+_fluxU0', 'O+_fluxU90', 'O+_fluxU180', 'O+_fluxPA', 'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 
                  'O+_fluxU_Dusk', 'e_fluxU', 'e_fluxU0', 'e_fluxU90', 'e_fluxU180', 'e_fluxPA']
    
        #Set up dictionary of colours for spacecrafts
        ColourDict = {}
        ColourDict.update({'C1': 'k'})
        ColourDict.update({'C2': 'r'})
        ColourDict.update({'C3': 'g'})
        ColourDict.update({'C4': 'b'})
    
        #If needed, use this to convert parameter 'codename' to real name on plot
        ParaName = {'Bmag': '|\mathbf{B}|', 'Bx': 'B_X', 'By': 'B_Y', 'Bz': 'B_Z', 'xgse': 'X_{GSE}', 'ygse': 'Y_{GSE}',
                    'zgse': 'Z_{GSE}', 'E': '|\mathbf{E}|', 'Ex': 'E_X', 'Ey': 'E_Y', 'Ez': 'E_Z', 'vi': '|\mathbf{v_i}|',
                    'vxi': 'v_{xi}', 'vyi': 'v_{yi}', 'vzi': 'v_{zi}', 'vH': '|\mathbf{v_H}|', 'vxH': 'v_{xH}',
                    'vyH': 'v_{yH}', 'vzH': 'v_{zH}', 'vHe': '|\mathbf{v_{He}}|', 'vxHe': 'v_{xHe}', 'vyHe': 'v_{yHe}', 'vzHe': 'v_{zHe}',
                    'vO': '|\mathbf{v_O}|', 'vxO': 'v_{xO}', 'vyO': 'v_{yO}', 'vzO': 'v_{zO}', 'Nion': 'N_i', 'Pion': 'P_i', 'Tion': 'T_i',
                    'TionBll': 'T_{i,\parallel}', 'TionBp': 'T_{i,\perp}', 'NH+': 'N_{H^+}', 'PH+': 'P_{H^+}', 'TH+': 'T_{H^+}',
                    'TH+Bll': 'T_{H^+,\parallel}', 'TH+Bp': 'T_{H^+,\perp}', 'NHe+': 'N_{He^+}', 'PHe+': 'P_{He^+}',
                    'THe+': 'T_{He^+}', 'THe+Bll': 'T_{He^+,\parallel}', 'THe+Bp': 'T_{He^+,\perp}', 'NO+': 'N_{O^+}',
                    'PO+': 'P_{O^+}', 'TO+': 'T_{O^+}', 'TO+Bll': 'T_{O^+,\parallel}', 'TO+Bp': 'T_{O^+,\perp}',
                    'Bvec': '\mathbf{B}', 'Evec': '\mathbf{E}', 'ion_vel': '\mathbf{v_i}', 'H+_vel': '\mathbf{v_H}',
                    'He+_vel': '\mathbf{v_He}', 'O+_vel': '\mathbf{v_O}', 'e_vel': '\mathbf{v_e}',
                    'ion_fluxU': 'Ions By-Energy Differential Energy Flux',
                    'ion_fluxU0': 'Ions By-Energy Differential Energy Flux (PA = 0^\circ)',
                    'ion_fluxU90': 'Ions By-Energy Differential Energy Flux (PA = 90^\circ)',
                    'ion_fluxU180': 'Ions By-Energy Differential Energy Flux (PA = 180^\circ)',
                    'ion_fluxPA': 'Ions By-PA Differential Energy Flux',
                    'ion_fluxU_Sun': 'Ions By-energy Differential Energy Flux (Sunward)',
                    'ion_fluxU_Dawn': 'Ions By-energy Differential Energy Flux (Dawnward)',
                    'ion_fluxU_Antisun': 'Ions By-energy Differential Energy Flux (Antisunward)',
                    'ion_fluxU_Dusk': 'Ions By-energy Differential Energy Flux (Duskward)',
                    'H+_fluxU': '$\mathrm{H^+}$ By-Energy Differential Energy Flux',
                    'H+_fluxU0': '$\mathrm{H^+}$ By-Energy Differential Energy Flux (PA = 0^\circ)',
                    'H+_fluxU90': '$\mathrm{H^+}$ By-Energy Differential Energy Flux (PA = 90^\circ)',
                    'H+_fluxU180': '$\mathrm{H^+}$ By-Energy Differential Energy Flux (PA = 180^\circ)',
                    'H+_fluxPA': '$\mathrm{H^+}$ By-PA Differential Energy Flux',
                    'H+_fluxU_Sun': '$\mathrm{H^+}$ By-energy Differential Energy Flux (Sunward)',
                    'H+_fluxU_Dawn': '$\mathrm{H^+}$ By-energy Differential Energy Flux (Dawnward)',
                    'H+_fluxU_Antisun': '$\mathrm{H^+}$ By-energy Differential Energy Flux (Antisunward)',
                    'H+_fluxU_Dusk': '$\mathrm{H^+}$ By-energy Differential Energy Flux (Duskward)',
                    'He+_fluxU': '$\mathrm{He^+}$ By-Energy Differential Energy Flux',
                    'He+_fluxU0': '$\mathrm{He^+}$ By-Energy Differential Energy Flux (PA = 0^\circ)',
                    'He+_fluxU90': '$\mathrm{He^+}$ By-Energy Differential Energy Flux (PA = 90^\circ)',
                    'He+_fluxU180': '$\mathrm{He^+}$ By-Energy Differential Energy Flux (PA = 180^\circ)',
                    'He+_fluxPA': '$\mathrm{He^+}$ By-PA Differential Energy Flux',
                    'He+_fluxU_Sun': '$\mathrm{He^+}$ By-energy Differential Energy Flux (Sunward)',
                    'He+_fluxU_Dawn': '$\mathrm{He^+}$ By-energy Differential Energy Flux (Dawnward)',
                    'He+_fluxU_Antisun': '$\mathrm{He^+}$ By-energy Differential Energy Flux (Antisunward)',
                    'He+_fluxU_Dusk': '$\mathrm{He^+}$ By-energy Differential Energy Flux (Duskward)',
                    'O+_fluxU': '$\mathrm{O^+}$ By-Energy Differential Energy Flux',
                    'O+_fluxU0': '$\mathrm{O^+}$ By-Energy Differential Energy Flux (PA = 0^\circ)',
                    'O+_fluxU90': '$\mathrm{O^+}$ By-Energy Differential Energy Flux (PA = 90^\circ)',
                    'O+_fluxU180': '$\mathrm{O^+}$ By-Energy Differential Energy Flux (PA = 180^\circ)',
                    'O+_fluxPA': '$\mathrm{O^+}$ By-PA Differential Energy Flux',
                    'O+_fluxU_Sun': '$\mathrm{O^+}$ By-energy Differential Energy Flux (Sunward)',
                    'O+_fluxU_Dawn': '$\mathrm{O^+}$ By-energy Differential Energy Flux (Dawnward)',
                    'O+_fluxU_Antisun': '$\mathrm{O^+}$ By-energy Differential Energy Flux (Antisunward)',
                    'O+_fluxU_Dusk': '$\mathrm{O^+}$ By-energy Differential Energy Flux (Duskward)',
                    'e_fluxU': 'Electrons By-Energy Differential Energy Flux',
                    'e_fluxU0': 'Electrons By-Energy Differential Energy Flux (PA = 0^\circ)',
                    'e_fluxU90': 'Electrons By-Energy Differential Energy Flux (PA = 90^\circ)',
                    'e_fluxU180': 'Electrons By-Energy Differential Energy Flux (PA = 180^\circ)',
                    'e_fluxPA': 'Electrons By-PA Differential Energy Flux'}
    
        #Remove panels with no data
        for key in PanelDict.keys():
            if np.size(PanelDict[key]) == 0:
                PanelDict.pop(key)

        #Number of subplot panels = number of keys in plotting dictionary
        no_subplot = len(list(PanelDict.keys()))
        paralist = list(PanelDict.keys())

        #Remove spacecraft potential as separate plot
        if 'SC_Pot' in paralist:
            no_subplot = no_subplot - 1
            
            #Record spacecraft potential data in CSV if required.
            if data_record == True:
                for sc in sclist:
                    comb_data = PanelDict['SC_Pot'][sc]
                    self.csv_record('SC_Pot', sc, comb_data)
                    print('{}, {} data recorded in csv'.format('SC_Pot',sc))

        if mode == 'light':
            ctick = 'k'
            cxline = 'r'
            cyline = 'k'
            czline = 'b'
            tchoice = False
        elif mode == 'dark':
            ctick = 'w'
            cxline = 'r'
            cyline = 'w'
            czline = 'cyan'
            tchoice = True
        else:
            raise Exception('Incorrect Colour Mode')
    
    
        #Set up figure and fontsizes
        fig, ax = plt.subplots(no_subplot, 1, figsize=(20, 5*no_subplot),squeeze=False)
        subtitle_fontsize = 32
        axlabel_fontsize = 24
        axtick_fontsize = 20
        cblabel_fontsize = 24
        cbtick_fontsize = 20
        legend_fontsize = 14

        #Plot each subplot
        for i in np.arange(no_subplot):

            para = list(PanelDict.keys())[i]
            
            if para in Multiple:

                if para == 'SC_Pot': #Remove spacecraft potential as separate plot
                    continue

                for sc in sclist:
                    comb_data = PanelDict[para][sc]

                    if np.size(comb_data) == 0:
                        continue
                    
                    dtime = comb_data[:, 0]

                    #Plot vector data in one panel
                    #Note: these variables exist only if veccomb is set to True
                    if para in ['Bvec', 'Evec', 'ion_vel', 'SCpos','e_vel']:  # vector
                        para_x = comb_data[:, 1]
                        para_y = comb_data[:, 2]
                        para_z = comb_data[:, 3]

                        #Remove fill values (9*10^30)
                        para_x[abs(para_x) > 10**10] = np.nan
                        para_y[abs(para_y) > 10**10] = np.nan
                        para_z[abs(para_z) > 10**10] = np.nan
                        
                        #Set label
                        if para == 'Bvec':
                            labelpara = 'B'
                        elif para == 'Evec':
                            labelpara = 'E'
                        elif para == 'ion_vel':
                            labelpara = 'Vi'
                        elif para == 'SCpos':
                            labelpara = 'R'
                        elif para == 'e_vel':
                            labelpara == 'Ve'

                        #Ignore spacecraft colour convention
                        ax[i,0].plot(dtime, para_x, label= 'X',
                                     color=cxline, ls='-', lw=2, alpha=1, zorder=2)
                        ax[i,0].plot(dtime, para_y, label= 'Y',
                                     color=cyline, ls='-', lw=2, alpha=1, zorder=2)
                        ax[i,0].plot(dtime, para_z, label= 'Z',
                                     color=czline, ls='-', lw=2, alpha=1, zorder=2)
                            
                        #Add line at zero if needed.
                        if para_x.min() < 0 and para_x.max() > 0:
                            ax[i,0].axhline(y=0,color=ctick,lw=1,alpha=1,zorder=1)
                        elif para_y.min() < 0 and para_y.max() > 0:
                            ax[i,0].axhline(y=0,color=ctick,lw=1,alpha=1,zorder=1)
                        elif para_z.min() < 0 and para_z.max() > 0:
                            ax[i,0].axhline(y=0,color=ctick,lw=1,alpha=1,zorder=1)

                    else:  #Plot scalar parameters.
                        para_data = comb_data[:, 1]

                        ax[i,0].plot(dtime, para_data, label=sc,
                                     color=ColourDict[sc], lw=2)
                        if para_data.min() < 0 and para_data.max() > 0:
                            ax[i,0].axhline(y=0,color=ctick,lw=1,alpha=1,zorder=1)
                            


                    
                #Highlight MP crossing times from Phan et al., 2003
                if MP_HL == True:
                    MP_start = dt.datetime(2002, 3, 18, 14, 54, 12)
                    MP_end = dt.datetime(2002, 3, 18, 14, 56, 10)
                    
                    #Check that start/end datetime not outside of input range
                    if MP_start < self.start_dt or MP_end > self.stop_dt: 
                        print(
                            'Magnetopause crossing interval not in plot range. No highlight.')
                    else:
                        if no_subplot > 1:
                            ax[i].axvspan(
                                MP_start, MP_end, facecolor='lightgray', edgecolor='None', alpha=0.8, zorder=0)
                        else:
                            ax.axvspan(MP_start, MP_end, facecolor='lightgray',
                                       edgecolor='None', alpha=0.8, zorder=0)
                
                #Highlight Walen test interval from Phan et al. 2003
                if Walen_HL == True:
                    Walen_start = dt.datetime(2002, 3, 18, 14, 57, 25)
                    Walen_end = dt.datetime(2002, 3, 18, 15, 3, 6)

                    if Walen_start < self.start_dt or Walen_end > self.stop_dt:
                        print('Walen Test interval not in plot range. No highlight.')
                    else:
                        if no_subplot > 1:
                            ax[i].axvline(Walen_start, color='purple',
                                          ls='--', lw=1.5, zorder=1)
                            ax[i].axvline(Walen_end, color='purple',
                                          ls='--', lw=1.5, zorder=1)
                        else:
                            ax.axvline(Walen_start, color='purple',
                                       ls='--', lw=1.5, zorder=1)
                            ax.axvline(Walen_end, color='purple',
                                       ls='--', lw=1.5, zorder=1)

                #Set up x-tick locations
                dtime_ticks = []
                dtime_ticklabels = []
                    
                #10 minutes major tick intervals for time period more than 20 minutes, 1 minute intervals for less.
                if (self.stop_dt-self.start_dt).seconds > 1200:
                    major_tick_interval = 600 #seconds
                else:
                    major_tick_interval = 60 #seconds
                    
                #In 1 or 10-minute intervals with decimals
                dtime_range = (self.stop_dt - self.start_dt).seconds/major_tick_interval
                m = 0
                
                while m <= dtime_range:
                    tick_dt = self.start_dt + m*dt.timedelta(seconds=600)
                    dtime_ticks.append(tick_dt)
                    tick_time_split = str(tick_dt).split(' ')[1].split(':')
                    tick_dt_str = ':'.join(
                        (tick_time_split[0], tick_time_split[1]))
                    dtime_ticklabels.append(tick_dt_str)
                    m = m+1

                #Set axes visiblity based on input.
                if xax_view == 0: #Set the x-axis to be visible for the bottommost subplot.
                    if i < no_subplot-1:
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
                        ax[i,0].set_xlim(self.start_dt, self.stop_dt)
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
                        ax[i,0].set_xlim(self.start_dt, self.stop_dt)
                        dt_str = self.start_dt.strftime(
                            '%Y-%m-%d')  # Set date string at the bottom of the figure.
                        ax[i,0].set_xlabel('{} (UT)'.format(
                            dt_str), fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif xax_view == 1: #Set the x-axis to be invisible for all subplots.
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
                    ax[i,0].set_xlim(self.start_dt, self.stop_dt)
                elif xax_view == 2: #Set the x-axis to be visible for all subplots.
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
                    ax[i,0].set_xlim(self.start_dt, self.stop_dt)
                    dt_str = self.start_dt.strftime(
                        '%Y-%m-%d')
                    ax[i,0].set_xlabel('{} (UT)'.format(
                        dt_str), fontsize=axlabel_fontsize, fontfamily='sans-serif')

                #Set y-axis: linear scale except B-field magnitude (log)
                #5 minor intervals between two major ticks
                if para not in ['Bmag']:
                    ax[i,0].yaxis.set_minor_locator(AutoMinorLocator(5))
                    ax[i,0].tick_params(axis='y', which='both',
                                      left=True, right=False)
                    ax[i,0].tick_params(axis='y', which='major',
                                      direction='inout', width=1, length=15, labelsize=axtick_fontsize)
                    ax[i,0].tick_params(axis='y', which='minor',
                                      direction='inout', width=1, length=5)

                #Set y-axis limits for different parameters.
                if para in ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'x', 'y', 'z']:

                    if para in ['Bx', 'By', 'Bz']:
                        ymin = -100
                        ymax = 100
                        ylabel = '{} (nT)'.format(para)
                    elif para in ['Ex', 'Ey', 'Ez']:
                        ymin = -100
                        ymax = 100
                        ylabel = '{} (mV/m)'.format(para)
                    elif para in ['x', 'y', 'z']:
                        ymin = -10**5
                        ymax = 10**5
                        ylabel = '{} (km)'.format(para)

                    ax[i,0].set_ylim([ymin, ymax])
                    ax[i,0].set_ylabel(ylabel, fontsize=axlabel_fontsize, fontfamily='sans-serif')
                    

                elif para == 'Bmag':
                    ymin = 1
                    ymax = 1000
                    ax[i,0].set_yscale('log')
                    ax[i,0].set_ylim([ymin, ymax])
                    ax[i,0].set_ylabel(
                        '$\mathrm{|\mathbf{B}|}$ (nT)', fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif para == 'ion_vel':
                    ax[i,0].set_ylabel(
                        '$\mathrm{v_i}$ (km/s)', fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif para == 'Nion':
                    ax[i,0].set_ylabel(
                        '$\mathrm{N_i}$ ($\mathrm{cm^{-3}}$)', fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif para == 'Pion':
                    ax[i,0].set_ylabel(
                        '$\mathrm{P_i}$ (nPa)', fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif para == 'Tion':
                    ax[i,0].set_ylabel(
                        '$\mathrm{T_i}$ (keV)', fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif para == 'Bvec':
                    ax[i,0].set_ylabel(
                        'B (nT)', fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif para == 'Evec':
                    ax[i,0].set_ylabel(
                        '$\mathrm{E_i}$ (mV/m)', fontsize=axlabel_fontsize, fontfamily='sans-serif')
                    

                #Spacecraft name text set up
                if len(sclist) > 1:
                    scstr = sclist[0]
                    for m in np.arange(1, len(sclist)):
                        scstr = scstr + ',{} '.format(sclist[m])
                else:
                    scstr = sclist[0]

                #Set in-plot labels
                #Subplot Label Letter
                ax[i,0].text(0.04, 0.85, '({})'.format(ALC[i+self.sp_ind0]), fontsize=subtitle_fontsize, fontfamily='sans-serif',
                             horizontalalignment='center', verticalalignment='center',
                             transform=ax[i,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                #Spacecraft name
                ax[i,0].text(0.96, 0.85, scstr, fontsize=subtitle_fontsize, fontfamily='sans-serif',
                             horizontalalignment='center', verticalalignment='center',
                             transform=ax[i,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                #Legend
                if len(sclist) > 1:
                    ax[i,0].legend(loc='lower right', prop={
                        'size': legend_fontsize,'family': 'sans-serif'})
                elif len(sclist) == 1 and para in ['Bvec', 'Evec', 'SCpos','ion_vel','e_vel']:
                    ax[i,0].legend(loc='lower right', prop={
                        'size': legend_fontsize,'family': 'sans-serif'})
                
                #Set colors
                ax[i,0].spines['bottom'].set_color(ctick) #Set axes color
                ax[i,0].spines['top'].set_color(ctick)
                ax[i,0].spines['left'].set_color(ctick)
                ax[i,0].spines['right'].set_color(ctick)
                ax[i,0].tick_params(axis='both', which='both', colors=ctick) #Set tick color
                ax[i,0].xaxis.label.set_color(ctick) #Set x-axis ticklabel color
                ax[i,0].yaxis.label.set_color(ctick) #Set y-axis ticklabel color
                
                #Save plot data to csv file.
                if data_record == True:
                    self.csv_record(para, sc, comb_data)
                    print('{}, {} data recorded in csv'.format(para, sc))


            #3D colour plot for flux data
            elif para in Single:
                
                #Only pick spacecraft with data.
                #Iterate with while loop until first available spacecraft is found
                no = 0
                switch = True

                while switch == True and no < len(sclist):
                    sc = sclist[no]
                    if np.size(PanelDict[para][sc]) == 0:
                        no = no+1
                    else:
                        switch = False

                if 'flux' in para: #Plotting function for flux data.
                    comb_data = PanelDict[para][sc]

                    if np.size(comb_data) == 0:  #Skip if no data.
                        continue
                    
                    #Remove NaN from datetime data.
                    dtime = comb_data[:, 0]
                    l1 = len(dtime)
                    dtime = np.array(
                        [d for d in dtime if type(d) is dt.datetime])
                    l2 = len(dtime)
                    
                    #Remove Nan from energy data and set up meshgrid.
                    if 'U' in para:
                        energy = comb_data[:, 1]
                        energy = np.array(
                            [e for e in energy if np.isnan(e) == False])
                        h, v = np.meshgrid(dtime, energy)
                    elif 'PA' in para:
                        pa = comb_data[:, 1]
                        pa = np.array([p for p in pa if np.isnan(p) == False])
                        h, v = np.meshgrid(dtime, pa)

                    #Remove nan from flux data.
                    if l1 > l2:  #If there has been reduction in dt length. Work only if nan is stacked at the end.
                        flux = comb_data[0:l2, 2:].transpose()
                    else:
                        flux = comb_data[:, 2:].transpose()

                    #Must convert flux array from object to float
                    fluxfloat = np.zeros(np.shape(flux))
                    for m in np.arange(len(flux[:, 0])):
                        for n in np.arange(len(flux[0, :])):
                            fluxfloat[m, n] = flux[m, n]
                    
                    #Set min/max for flux (can depend on datetime range used)
                    if para in ['ion_fluxU']:  # set min/max for count
                        flux_min = 10**4  
                        flux_max = 10**9
                        vert_min = 0.01
                        vert_max = energy.max()
                        vert_label = 'Energy (keV)'
                        cblabel = para.split('_')[0].title(
                        ) + ' DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                    elif para in ['H+_fluxU', 'He+_fluxU', 'O+_fluxU']:
                        flux_min = 10**3  
                        flux_max = 10**8
                        vert_min = 0.01
                        vert_max = energy.max()
                        vert_label = 'Energy (keV)'
                        cblabel = para.split('_')[0].title(
                        ) + ' DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                    elif para in ['ion_fluxU0', 'ion_fluxU90', 'ion_fluxU180']:
                        flux_min = 10**6
                        flux_max = 10**9
                        vert_min = 5
                        vert_max = energy.max()
                        vert_label = 'Energy (keV)'
                        cblabel = para.split('_')[0].title(
                        ) + ' DEF (PA = {})'.format(para.split('U')[-1]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                    elif para in ['H+_fluxU0', 'H+_fluxU90', 'H+_fluxU180',
                                  'He+_fluxU0', 'He+_fluxU90', 'He+_fluxU180',
                                  'O+_fluxU0', 'O+_fluxU90', 'O+_fluxU180']:
                        flux_min = 10**3  
                        flux_max = 10**6 
                        vert_min = 0.01
                        vert_max = energy.max()
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
                        vert_max = energy.max()
                        vert_label = 'Energy (keV)'
                        cblabel = para.split('_')[0].title(
                        ) + ' DEF ({}ward)'.format(para.split('_')[2]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                    elif para in ['H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk',
                                  'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk',
                                  'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk']:
                        flux_min = 10**3
                        flux_max = 10**5
                        vert_min = 0.05
                        vert_max = energy.max()
                        vert_label = 'Energy (keV)'
                        cblabel = para.split('_')[0].title(
                        ) + ' DEF ({}ward)'.format(para.split('_')[2]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                    elif para in ['e_fluxU']:
                        flux_min = 5*10**6
                        flux_max = 5*10**9
                        vert_min = 5/1000
                        vert_max = energy.max()
                        vert_label = 'Energy (keV)'
                        cblabel = 'Electron DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                    elif para in ['e_fluxU0', 'e_fluxU90', 'e_fluxU180']:
                        flux_min = 5*10**6
                        flux_max = 3*10**9
                        vert_min = 5/1000
                        vert_max = energy.max()
                        vert_label = 'Energy (keV)'
                        cblabel = 'Electron DEF (PA = {})'.format(para.split('U')[-1]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                    elif para == 'e_fluxPA':
                        flux_min = 10**6  
                        flux_max = 10**8
                        vert_min = 0
                        vert_max = 180
                        vert_label = 'Pitch Angle ($^\circ$)'
                        cblabel = 'Electron DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'

                    #Set color map
                    cmap_name = 'nipy_spectral'
                    
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
                        cm = colors.LinearSegmentedColormap.from_list(cmap_name2,clist2,n_bins) #Create new colormap and set it to this name
                    else:
                        cm = plt.get_cmap(cmap_name) #In dark mode, keep using nipy_spectral.
                    
                    #Plot fkux data.
                    im = ax[i,0].pcolormesh(h, v, fluxfloat, norm=colors.LogNorm(
                        vmin=flux_min, vmax=flux_max), cmap=cm, edgecolors='none')
                    

                    #If spacecraft potential has been calculated, applied to all U0, U90, U180 e flux plots
                    if 'SC_Pot' in paralist and 'e_fluxU' in para:
                        photoe_dtime = PanelDict['SC_Pot'][sc][:, 0]
                        photoe_energy = -1 * \
                            PanelDict['SC_Pot'][sc][:, 1]  #Convert to e
                        photoe_energy = photoe_energy/1000 #Convert to keV
                            #Keep line black for clarity.
                        ax[i,0].plot(photoe_dtime, photoe_energy,
                                     c='k',ls='solid', linewidth=2)
                        #Add line for 30% increase used as energy threshold.
                        ax[i,0].plot(photoe_dtime, photoe_energy *
                                     1.30, c='k',ls=(5,(10,5)), linewidth=2) #long dash with offset


                    #Set x-axis properties. See notes in multiple parameter section
                    dtime_ticks = []  
                    dtime_ticklabels = []
                    
                    if (self.stop_dt-self.start_dt).seconds > 1200:
                        major_tick_interval = 600 
                    else:
                        major_tick_interval = 60 
                    
                    dtime_range = (self.stop_dt - self.start_dt).seconds/major_tick_interval
                    m = 0
                    while m <= dtime_range:
                        tick_dt = self.start_dt + m*dt.timedelta(seconds=major_tick_interval)
                        dtime_ticks.append(tick_dt)
                        tick_time_split = str(tick_dt).split(' ')[1].split(':')
                        tick_dt_str = ':'.join(
                            (tick_time_split[0], tick_time_split[1]))
                        dtime_ticklabels.append(tick_dt_str)
                        m = m+1

                    if xax_view == 0:
                        if i < no_subplot-1:
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
                            ax[i,0].set_xlim(self.start_dt, self.stop_dt)
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
                            ax[i,0].set_xlim(self.start_dt, self.stop_dt)
                            dt_str = self.start_dt.strftime(
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
                        ax[i,0].set_xlim(self.start_dt, self.stop_dt)
                    elif xax_view == 2:
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
                        ax[i,0].set_xlim(self.start_dt, self.stop_dt)
                        dt_str = self.start_dt.strftime(
                            '%Y-%m-%d')
                        ax[i,0].set_xlabel('{} (UT)'.format(
                            dt_str), fontsize=axlabel_fontsize, fontfamily='sans-serif')

                    #Set y-axis properties, differing between flux sorted by energy (U) and by pitch angle (PA)
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
                
                
                    #Set subplot label letter.
                    ax[i,0].text(0.04, 0.85, '({})'.format(ALC[i+self.sp_ind0]), fontsize=subtitle_fontsize, fontfamily='sans-serif',
                                 horizontalalignment='center', verticalalignment='center',
                                 transform=ax[i,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                    #Set spacecraft label.
                    ax[i,0].text(0.96, 0.85, sc, fontsize=subtitle_fontsize, fontfamily='sans-serif',
                                 horizontalalignment='center', verticalalignment='center',
                                 transform=ax[i,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                    
                    #Set subplot title if desired.
                    #ax[i,0].title.set_text('({}) {} {}'.format(ALC[i+self.sp_ind0],ParaName[para],sc))
                    #ax[i,0].title.set_fontsize(subtitle_fontsize)
                        
                    #Set subplot color properties.
                    ax[i,0].spines['bottom'].set_color(ctick)
                    ax[i,0].spines['top'].set_color(ctick)
                    ax[i,0].spines['left'].set_color(ctick)
                    ax[i,0].spines['right'].set_color(ctick)
                    ax[i,0].tick_params(axis='both',which='both',colors=ctick) 
                    ax[i,0].xaxis.label.set_color(ctick) 
                    ax[i,0].yaxis.label.set_color(ctick)
                    cb.ax.axes.tick_params(axis='y',which='both',colors=ctick) 
                    cb.ax.axes.yaxis.label.set_color(ctick) 
                    
                    #Record plot data in CSV file if needed.   
                    if data_record == True:
                        self.csv_record(para, sc, comb_data)
                        print('{}, {} data recorded in csv'.format(para, sc))
        
        
        #Save figure
        targetfolder = 'ClusterPanels/{}/{}/{}/'.format(self.start_dt.year,self.start_dt.month,self.start_dt.day)
        os.makedirs(targetfolder,exist_ok=True)
        start_str = dt.datetime.strftime(self.start_dt, '%Y%m%d%H%M%S')
        stop_str = dt.datetime.strftime(self.stop_dt, '%Y%m%d%H%M%S')
        savename = targetfolder + '{}_{}_'.format(start_str,stop_str)

        for sc in sclist:  # add sc to savename
            savename = savename + '{}'.format(sc)

        for para in paralist:  # add para to savename
            savename = savename + '_{}'.format(para)

        savename = savename + '.{}'.format(saveformat)  # add file extension
        plt.subplots_adjust(left=0.08, bottom=0.05, right=0.90, top=0.97,
                            wspace=0.0, hspace=0.08)  # 0.05 = 5% from left figure width/height
        #Fig 3 ([0.08,0.05,0.90,0.97,0.0,0.08])
        #Fig 4 ([0.08,0.05,0.90,0.97,0.0,0.08])
        
        plt.savefig(savename, dpi=200, format=saveformat, transparent=tchoice)
        plt.close()
        
    def PanelCompArrange(self, paralist, sclist, DataDict):
        #Create plotting dictionary in which all four spacecraft flux data is compared.
        
        PanelCompDict = {}

        #Only 'single' parameters
        Single = ['ion_fluxU', 'ion_fluxU0', 'ion_fluxU90', 'ion_fluxU180', 'ion_fluxPA', 
                  'ion_fluxU_Sun', 'ion_fluxU_Dawn', 'ion_fluxU_Antisun', 'ion_fluxU_Dusk', 
                  'H+_fluxU', 'H+_fluxU0', 'H+_fluxU90', 'H+_fluxU180', 'H+_fluxPA',
                  'H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk', 
                  'He+_fluxU', 'He+_fluxU0', 'He+_fluxU90', 'He+_fluxU180', 'He+_fluxPA', 
                  'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk', 
                  'O+_fluxU', 'O+_fluxU0', 'O+_fluxU90', 'O+_fluxU180', 'O+_fluxPA', 
                  'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk', 
                  'e_fluxU', 'e_fluxU0', 'e_fluxU90', 'e_fluxU180', 'e_fluxPA']
        
        Para3D = ['ion_fluxU_Sun', 'ion_fluxU_Dawn', 'ion_fluxU_Antisun', 'ion_fluxU_Dusk', 
                    'H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk', 
                    'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk', 
                    'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk']

        #See PanelArrange function for comments. Same principles are utilised
        #EXCEPT no while loop. All spacecrafts are plotted. If there is no data, plot shouldl be blank.

        for para in paralist:

            if para not in Single:  #Skip inapplicable parameter
                continue

            PanelCompDict.update({para: {}})
            for sc in sclist:
                PanelCompDict[para].update({sc: []})

            for key in DataDict.keys():
                for sc in sclist:
                    convert_in = {sc: [para]}
                    convert_out = self.VarSorter(
                        convert_in, self.instrlist, prior_sorted=True)

                    ID = convert_out[sc][0].split(
                        '__')[1]

                    if ID == key:
                        print('Key {} ID Matched'.format(ID))

                        dt_arr = np.empty([1])
                        energy_arr = np.empty([1])
                        pa_arr = np.empty([1])

                        if 'CIS' in ID:
                            if 'U' in para:
                                #16 = pa, 31 = energy
                                data_arr = np.empty([1, 31])
                            elif 'PA' in para:
                                data_arr = np.empty([1, 16])
                        elif 'PEA' in ID:
                            if 'U' in para:
                                #12 = pa, 44 = energy
                                data_arr = np.empty([1, 44])
                            elif 'PA' in para:
                                data_arr = np.empty([1, 12])

                        for period in DataDict[ID].keys():

                            dtime = DataDict[ID][period]['datetime']
                            if isinstance(dtime[0], dt.datetime) == False:
                                print(
                                    'Append Stop Due to Datetime Issue for ID {} and Period {}'.format(ID, period))
                                continue

                            energy = DataDict[ID][period]['energy']
                            if 'PEA' in ID:
                                if len(np.shape(energy)) != 2:  
                                    continue
                                energy = energy[0, :]/1000
                            else:
                                if 'fluxU' in para:
                                    energy = energy/1000  
                                else:
                                    energy = energy  

                            if para in Para3D: #3D distributions
                                theta = DataDict[ID][period]['theta']
                                phi = DataDict[ID][period]['phi']
                            else: #PAD distributions
                                pa = DataDict[ID][period]['pa']

                            pdata_init = np.copy(DataDict[ID][period]['diffflux'])
                            if para in Para3D:
                                if len(np.shape(pdata_init)) != 4: 
                                    continue
                            else:
                                if len(np.shape(pdata_init)) != 3: 
                                    continue

                            #For CIS, convert particle flux to energy flux
                            if 'CIS' in ID:
                                if para in Para3D:
                                    for i in np.arange(len(dtime)):
                                        for j in np.arange(len(theta)):
                                            for k in np.arange(len(phi)):
                                                pdata_init[i,j,k,:] = pdata_init[i,j,k,:]*energy
                                    
                                else:
                                    for i in np.arange(len(dtime)):
                                        for j in np.arange(len(pa)):
                                            if 'fluxU' in para:
                                                pdata_init[i,j,:] = pdata_init[i, j,:]*energy
                                            else:
                                                pdata_init[i,j,:] = pdata_init[i,j,:]*energy/1000

                            if 'fluxU' in para:  #Flux by energy from pitch angle distributions
                                pdata = np.zeros((len(dtime), len(energy)))

                                for i in np.arange(len(dtime)):
                                    for j in np.arange(len(energy)):
                                        if 'U0' in para:
                                            pdata[i, j] = pdata_init[i, 0, j]
                                        elif 'U90' in para:
                                            if 'CIS' in ID:
                                                pdata_pa = pdata_init[i, 7:9, j]
                                            elif 'PEA' in ID:
                                                pdata_pa = pdata_init[i, 5:7, j]
                                            pdata[i, j] = np.sum(
                                                pdata_pa)/len(pdata_pa)
                                        elif 'U180' in para:
                                            pdata[i, j] = pdata_init[i, -1, j]
                                        else:
                                            pdata_pa = pdata_init[i, :, j]
                                            pdata[i, j] = np.sum(
                                                pdata_pa)/len(pdata_pa)

                                dt_arr = np.concatenate((dt_arr, dtime))
                                energy_arr = np.concatenate(
                                    (energy_arr, energy))
                                # should get dt x energy
                                data_arr = np.concatenate(
                                    (data_arr, pdata), axis=0)

                            elif 'fluxPA' in para: #Flux by pitch angle from pitch angle distributions
                                pdata = np.zeros((len(dtime), len(pa)))

                                for i in np.arange(len(dtime)):
                                    for j in np.arange(len(pa)):
                                        pdata_U = pdata_init[i, j, :]

                                        if 'SC_Pot' in paralist and 'e_flux' in para:  #Remove photo-electrons

                                            PotID = '{}_CP_EFW_L3_P'.format(
                                                sc.upper())
                                            Pot_Period = None
                                            for c in DataDict[PotID].keys():
                                                PP_startstr = c.split(
                                                    '___')[0]  
                                                PP_start = dt.datetime.strptime(
                                                    PP_startstr, '%Y-%m-%dT%H-%M-%S')
                                                PP_endstr = c.split('___')[1]
                                                PP_end = dt.datetime.strptime(
                                                    PP_endstr, '%Y-%m-%dT%H-%M-%S')

                                                if dtime[i] > PP_start and dtime[i] < PP_end:
                                                    Pot_Period = c

                                            if Pot_Period == None:
                                                print(
                                                    'Photo-Electron Energy Threshold not Found')
                                                pdata[i, j] = np.sum(
                                                    pdata_U)/len(pdata_U)
                                                continue

                                            Pot_dtime = DataDict[PotID][Pot_Period]['datetime']
                                            diff_arr = abs(dtime[i]-Pot_dtime)
                                            match_ind = np.where(
                                                diff_arr == min(diff_arr))
                                            match_Pot = DataDict[PotID][Pot_Period]['SC_Pot'][match_ind]

                                            photoe = -1*match_Pot/1000 
                                            photoe = photoe*1.30  
                                            if len(photoe) > 1:
                                                photoe = max(photoe)
                                            pdata_Unew = pdata_U[energy > photoe]
                                            if len(pdata_Unew) != 0:
                                                pdata[i, j] = np.sum(
                                                    pdata_Unew)/len(pdata_U)
                                            elif len(pdata_Unew) == 0:
                                                pdata[i, j] = np.sum(
                                                    pdata_U)/len(pdata_U)
                                        else:
                                            pdata[i, j] = np.sum(
                                                pdata_U)/len(pdata_U)

                                dt_arr = np.concatenate((dt_arr, dtime))
                                pa_arr = np.concatenate((pa_arr, pa))
                                data_arr = np.concatenate(
                                    (data_arr, pdata), axis=0)  

                            elif para in Para3D: #Flux by energy from 3D particle distributions
                                pdata = np.zeros((len(dtime), len(energy)))
                                for i in np.arange(len(dtime)):
                                    for j in np.arange(len(energy)):
                                        
                                        #Group azimuthal angles (accounted for direction of arrival)
                                        if 'fluxU_Sun' in para:
                                            pdata_th = pdata_init[i,:,2:10,j]
                                        elif 'fluxU_Dawn' in para:
                                            pdata_th = pdata_init[i,:,6:14,j]
                                        elif 'fluxU_Antisun' in para:
                                            pdata_thA = pdata_init[i,:,10:,j] 
                                            pdata_thB = pdata_init[i,:,0:2,j] 
                                            pdata_th = np.concatenate((pdata_thA,pdata_thB),axis=1)
                                        elif 'fluxU_Dusk' in para:
                                            pdata_thA = pdata_init[i,:,-2:,j] 
                                            pdata_thB = pdata_init[i,:,0:6,j]
                                            pdata_th = np.concatenate((pdata_thA,pdata_thB),axis=1)
                                        
                                        #Average over all azimuthal angles
                                        pdata_th = np.sum(pdata_th,axis=1)/len(pdata_th[0,:])
                                        
                                        #Then average over all elevation angles
                                        pdata[i,j] = np.sum(pdata_th)/len(pdata_th)
                                        
            
                        #Remove initialising first element
                        dt_arr = np.delete(dt_arr, (0), axis=0)
                        energy_arr = np.delete(energy_arr, (0), axis=0)
                        pa_arr = np.delete(pa_arr, (0), axis=0)
                        data_arr = np.delete(data_arr, (0), axis=0)

                        #Fill in NaN elements ahead of array combination
                        if 'fluxU' in para:
                            if len(dt_arr) >= len(energy_arr):
                                nan_extra = np.nan * \
                                    np.zeros((len(dt_arr)-len(energy_arr)))
                                energy_arr = np.concatenate(
                                    (energy_arr, nan_extra))
                            else:
                                nan_extra = np.nan * \
                                    np.zeros((len(energy_arr)-len(dt_arr)))
                                nan_extra2D = np.nan * \
                                    np.zeros(
                                        (len(energy_arr)-len(dt_arr), len(data_arr[0, :])))

                                dt_arr = np.concatenate((dt_arr, nan_extra))
                                data_arr = np.concatenate(
                                    (data_arr, nan_extra2D), axis=0)


                            comb1 = np.stack(
                                (dt_arr, energy_arr)).transpose()
                            comb2 = np.concatenate((comb1, data_arr), axis=1)

                        elif 'fluxPA' in para:

                            if len(dt_arr) >= len(pa_arr):
                                nan_extra = np.nan * \
                                    np.zeros((len(dt_arr)-len(pa_arr)))
                                pa_arr = np.concatenate((pa_arr, nan_extra))
                            else:
                                nan_extra = np.nan * \
                                    np.zeros((len(pa_arr)-len(dt_arr)))
                                nan_extra2D = np.nan * \
                                    np.zeros(
                                        (len(pa_arr)-len(dt_arr), len(data_arr[0, :])))

                                dt_arr = np.concatenate((dt_arr, nan_extra))
                                data_arr = np.concatenate(
                                    (data_arr, nan_extra2D), axis=0)

                            comb1 = np.stack(
                                (dt_arr, pa_arr)).transpose()  # make 2D
                            comb2 = np.concatenate((comb1, data_arr), axis=1)

                        #if there is no data, fill array with NaNs
                        if np.size(data_arr) == 0:
                            comb_nan = np.zeros(np.shape(comb2))
                            comb_nan[:] = np.nan
                            PanelCompDict[para][sc] = comb_nan
                        else:
                            PanelCompDict[para][sc] = comb2

        return PanelCompDict

    def PanelCompArrange(self, PanelCompDict, sclist, xax_view, mode, saveformat,data_record):
        #Create plotting dictionary in which all four spacecraft flux data is compared.
        
        PanelCompDict = {}

        #Only 'single' parameters
        Single = ['ion_fluxU', 'ion_fluxU0', 'ion_fluxU90', 'ion_fluxU180', 'ion_fluxPA', 
                  'ion_fluxU_Sun', 'ion_fluxU_Dawn', 'ion_fluxU_Antisun', 'ion_fluxU_Dusk', 
                  'H+_fluxU', 'H+_fluxU0', 'H+_fluxU90', 'H+_fluxU180', 'H+_fluxPA',
                  'H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk', 
                  'He+_fluxU', 'He+_fluxU0', 'He+_fluxU90', 'He+_fluxU180', 'He+_fluxPA', 
                  'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk', 
                  'O+_fluxU', 'O+_fluxU0', 'O+_fluxU90', 'O+_fluxU180', 'O+_fluxPA', 
                  'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk', 
                  'e_fluxU', 'e_fluxU0', 'e_fluxU90', 'e_fluxU180', 'e_fluxPA']
        
        Para3D = ['ion_fluxU_Sun', 'ion_fluxU_Dawn', 'ion_fluxU_Antisun', 'ion_fluxU_Dusk', 
                    'H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk', 
                    'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk', 
                    'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk']
        
        #See PanelArrange function for comments. Same principles are utilised
        #EXCEPT no while loop. All spacecrafts are plotted as separate figures. If there is no data, plot should be blank.

        subtitle_fontsize = 32
        axlabel_fontsize = 24
        axtick_fontsize = 20
        cblabel_fontsize = 24
        cbtick_fontsize = 20

        paralist = list(PanelCompDict.keys())

        for para in paralist:

            to_pop = []
            for sc in PanelCompDict[para].keys():
                if np.size(PanelCompDict[para][sc]) == 0:
                    to_pop.append(sc)

            for tp in to_pop:
                PanelCompDict[para].pop(tp) 

            no_subplot = len(list(PanelCompDict[para].keys()))
            sclist = list(PanelCompDict[para].keys())
            fig, ax = plt.subplots(no_subplot, 1, figsize=(20, 5*no_subplot),squeeze=False)

            for i in np.arange(no_subplot):
                sc = sclist[i]

                comb_data = PanelCompDict[para][sc]

                dtime = comb_data[:, 0]
                l1 = len(dtime)
                dtime = np.array([d for d in dtime if type(d) is dt.datetime])
                l2 = len(dtime)

                if 'U' in para:
                    energy = comb_data[:, 1]
                    energy = np.array(
                        [e for e in energy if np.isnan(e) == False])
                    h, v = np.meshgrid(dtime, energy)
                elif 'PA' in para:
                    pa = comb_data[:, 1]
                    pa = np.array([p for p in pa if np.isnan(p) == False])
                    h, v = np.meshgrid(dtime, pa)

                if l1 > l2:
                    flux = comb_data[0:l2, 2:].transpose()
                else:
                    flux = comb_data[:, 2:].transpose()

                fluxfloat = np.zeros(np.shape(flux))
                for m in np.arange(len(flux[:, 0])):
                    for n in np.arange(len(flux[0, :])):
                        fluxfloat[m, n] = flux[m, n]

                if para in ['ion_fluxU']:
                    flux_min = 10**4  
                    flux_max = 10**9  
                    vert_min = 0.01
                    vert_max = energy.max()
                    vert_label = 'Energy (keV)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['H+_fluxU', 'He+_fluxU', 'O+_fluxU']:
                    flux_min = 10**3  
                    flux_max = 10**6  
                    vert_min = 0.01
                    vert_max = energy.max()
                    vert_label = 'Energy (keV)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['ion_fluxU0', 'ion_fluxU90', 'ion_fluxU180']:
                    flux_min = 10**6
                    flux_max = 10**9
                    vert_min = 5
                    vert_max = energy.max()
                    vert_label = 'Energy (keV)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF (PA = {})'.format(para.split('U')[-1]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['H+_fluxU0', 'H+_fluxU90', 'H+_fluxU180',
                              'He+_fluxU0', 'He+_fluxU90', 'He+_fluxU180',
                              'O+_fluxU0', 'O+_fluxU90', 'O+_fluxU180']:
                    flux_min = 10**3  
                    flux_max = 10**7  
                    vert_min = 0.01
                    vert_max = energy.max()
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
                    vert_max = energy.max()
                    vert_label = 'Energy (keV)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF ({}ward)'.format(para.split('_')[2]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['H+_fluxU_Sun', 'H+_fluxU_Dawn', 'H+_fluxU_Antisun', 'H+_fluxU_Dusk',
                              'He+_fluxU_Sun', 'He+_fluxU_Dawn', 'He+_fluxU_Antisun', 'He+_fluxU_Dusk',
                              'O+_fluxU_Sun', 'O+_fluxU_Dawn', 'O+_fluxU_Antisun', 'O+_fluxU_Dusk']:
                    flux_min = 10**3
                    flux_max = 10**5
                    vert_min = 0.05
                    vert_max = energy.max()
                    vert_label = 'Energy (keV)'
                    cblabel = para.split('_')[0].title(
                    ) + ' DEF ({}ward)'.format(para.split('_')[2]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'    
                elif para in ['e_fluxU']:
                    flux_min = 5*10**6
                    flux_max = 5*10**9
                    vert_min = 5/1000
                    vert_max = energy.max()
                    vert_label = 'Energy (keV)'
                    cblabel = 'Electron DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para in ['e_fluxU0', 'e_fluxU90', 'e_fluxU180']:
                    flux_min = 5*10**6
                    flux_max = 5*10**9
                    vert_min = 5/1000
                    vert_max = energy.max()
                    vert_label = 'Energy (keV)'
                    cblabel = 'Electron DEF (PA = {})'.format(para.split('U')[-1]) + ' \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                elif para == 'e_fluxPA':
                    flux_min = 10**6 
                    flux_max = 10**9 
                    vert_min = 0
                    vert_max = 180
                    vert_label = 'Pitch Angle ($^\circ$)'
                    cblabel = 'Electron DEF \n $\mathrm{keV/(cm^2 \cdot s \cdot sr \cdot keV)}$'
                    
                if mode == 'light':
                    ctick = 'k'
                    tchoice = False
                elif mode == 'dark':
                    ctick = 'w'
                    tchoice = True
                else:
                    raise Exception('Incorrect Colour Mode')
                    
                cmap_name = 'nipy_spectral'
                if mode == 'light':
                    cm = plt.get_cmap(cmap_name,10).copy() 
                    clist = [] #colorlist
                    for a in range(10):
                        clist.append(cm(a)) 
                    clist2 = clist.copy()
                    clist2[0] = clist[-1] 
                    clist2[-1] = clist[0]
                    n_bins = 256 
                    cmap_name2 = 'nipy_swapped' 
                    cm = colors.LinearSegmentedColormap.from_list(cmap_name2,clist2,n_bins)
                else:
                    cm = plt.get_cmap(cmap_name) 

                im = ax[i,0].pcolormesh(h, v, fluxfloat, norm=colors.LogNorm(
                    vmin=flux_min, vmax=flux_max), cmap=cm, edgecolors='none')

                if 'SC_Pot' in self.paralist and 'e_fluxU' in para:
                    photoe_dtime = self.PanelDict['SC_Pot'][sc][:, 0]
                    photoe_energy = -1 * \
                        self.PanelDict['SC_Pot'][sc][:, 1] 
                    photoe_energy = photoe_energy/1000
                    ax[i,0].plot(photoe_dtime, photoe_energy,
                                 c='k',ls='solid', linewidth=2)
                    ax[i,0].plot(photoe_dtime, photoe_energy *
                                 1.30, c='k',ls=(5,(10,3)), linewidth=2)


                dtime_ticks = []
                dtime_ticklabels = []
                
                if (self.stop_dt-self.start_dt).seconds > 1200:
                    major_tick_interval = 600 
                else:
                    major_tick_interval = 60 

                dtime_range = (self.stop_dt - self.start_dt).seconds/major_tick_interval
                m = 0
                while m <= dtime_range:
                    tick_dt = self.start_dt + m*dt.timedelta(seconds=600)
                    dtime_ticks.append(tick_dt)
                    tick_time_split = str(tick_dt).split(' ')[1].split(':')
                    tick_dt_str = ':'.join(
                        (tick_time_split[0], tick_time_split[1]))
                    dtime_ticklabels.append(tick_dt_str)
                    m = m+1

                if xax_view == 0:
                    if i < no_subplot-1:
                        ax[i,0].minorticks_on()
                        ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(10))
                        ax[i,0].tick_params(
                            axis='x', which='both', bottom=True, top=False, labelbottom=False)
                        ax[i,0].tick_params(
                            axis='x', which='major', direction='inout', width=1, length=15)
                        ax[i,0].tick_params(
                            axis='x', which='minor', direction='inout', width=1, length=5)
                        ax[i,0].axes.get_xaxis().set_ticks(
                            dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize) 
                        ax[i,0].set_xlim(self.start_dt, self.stop_dt)
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
                        ax[i,0].set_xlim(self.start_dt, self.stop_dt)
                        dt_str = self.start_dt.strftime(
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
                        ax[i,0].tick_params(
                            axis='x', which='minor', direction='inout', width=1, length=5)
                        ax[i,0].axes.get_xaxis().set_ticks(
                            dtime_ticks, labels=dtime_ticklabels,fontsize=axtick_fontsize) 
                        ax[i,0].set_xlim(self.start_dt, self.stop_dt)
                elif xax_view == 2:
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
                    ax[i,0].set_xlim(self.start_dt, self.stop_dt)
                    dt_str = self.start_dt.strftime(
                        '%Y-%m-%d')  # set date string
                    ax[i,0].set_xlabel('{} (UT)'.format(
                        dt_str), fontsize=axlabel_fontsize, fontfamily='sans-serif')

                if 'U' in para:
                    ax[i,0].set_yscale('log')
                    ax[i,0].tick_params(
                        axis='y', which='major', direction='inout', width=1, length=15)
                    ax[i,0].tick_params(
                        axis='y', which='minor', direction='inout', width=1, length=5)
                    ax[i,0].set_ylim(vert_min, vert_max)
                    ax[i,0].set_ylabel(vert_label, fontsize=axlabel_fontsize, fontfamily='sans-serif')
                elif 'PA' in para:
                    ax[i,0].yaxis.set_minor_locator(AutoMinorLocator(5))
                    ax[i,0].tick_params(
                        axis='y', which='major', direction='inout', width=1, length=15)
                    ax[i,0].tick_params(
                        axis='y', which='minor', direction='inout', width=1, length=5)
                    ax[i,0].set_ylim(vert_min, vert_max)
                    ax[i,0].set_ylabel(vert_label, fontsize=axlabel_fontsize, fontfamily='sans-serif')


                axins = inset_axes(
                    ax[i,0], width="1%", height="100%", loc='right', borderpad=-2)
                cb = fig.colorbar(
                    im, cax=axins, orientation='vertical', pad=0.01)
                cb.set_label(label=cblabel, size=cblabel_fontsize, fontfamily='sans-serif')
                cb.ax.tick_params(labelsize=cbtick_fontsize)
                ax[i,0].tick_params(axis='both', which='major',
                                    labelsize=axtick_fontsize)
                ax[i,0].text(0.04, 0.85, '({})'.format(ALC[i+self.sp_ind0]), fontsize=subtitle_fontsize, fontfamily='sans-serif',
                             horizontalalignment='center', verticalalignment='center',
                             transform=ax[i,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                ax[i,0].text(0.96, 0.85, sc, fontsize=subtitle_fontsize, fontfamily='sans-serif',
                             horizontalalignment='center', verticalalignment='center',
                             transform=ax[i,0].transAxes, bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))
                
                #ax[i,0].title.set_text('({}) {} {}'.format(ALC[i+self.sp_ind0],ParaName[para],sc))
                #ax[i,0].title.set_fontsize(subtitle_fontsize)
                           
                ax[i,0].spines['bottom'].set_color(ctick)
                ax[i,0].spines['top'].set_color(ctick)
                ax[i,0].spines['left'].set_color(ctick)
                ax[i,0].spines['right'].set_color(ctick)
                ax[i,0].tick_params(axis='both',which='both',colors=ctick) 
                ax[i,0].xaxis.label.set_color(ctick)
                ax[i,0].yaxis.label.set_color(ctick)
                cb.ax.axes.tick_params(axis='y',which='both',colors=ctick) 
                cb.ax.axes.yaxis.label.set_color(ctick) 
                
                if data_record == True:
                    self.csv_record(para, sc, comb_data)
                    print('{}, {} data recorded in csv'.format(para, sc))

            #Save Figure
            targetfolder = 'ClusterPanels/{}/{}/{}/'.format(self.start_dt.year,self.start_dt.month,self.start_dt.day)
            os.makedirs(targetfolder,exist_ok=True)
            start_str = dt.datetime.strftime(self.start_dt, '%Y%m%d%H%M%S')
            stop_str = dt.datetime.strftime(self.stop_dt, '%Y%m%d%H%M%S')
            savename = targetfolder + '{}_{}_'.format(start_str,stop_str)

            for sc in sclist:  # add sc to savename
                savename = savename + '{}'.format(sc)

            for para in paralist:  # add para to savename
                savename = savename + '_{}'.format(para)

            savename = savename + 'SCComp.{}'.format(saveformat)  # add file extension
            plt.subplots_adjust(left=0.08, bottom=0.05, right=0.90, top=0.97,
                                wspace=0.0, hspace=0.08)  # 0.05 = 5% from left figure width/height

            #Fig 3 ([0.08,0.05,0.92,0.97,0.0,0.08])
            #Fig 4 ([0.08,0.03,0.92,0.97,0.0,0.08])
            plt.savefig(savename, dpi=200, format=saveformat, transparent=tchoice)
            plt.close()
            
    def GSM_from_GSE(self,dt_arr,data_arr):
        #Convert GSE data array to GSM
        #To be called by PanelArrange, after period-stacking but before combining with datetime/energy.
        #Take datetime arrays because conversion is datetime dependent.
        
        #Make sure that data_arr has datetime along rows.
        #Create SpacePy GSE object then convert to GSM.
        GSE_obj = coord.Coords(data_arr,'GSE','car')
        GSE_obj.ticks = Ticktock(dt_arr,'UTC') #datetime on CSA is a UTC
        GSM_obj = GSE_obj.convert('GSM','car')
        
        return GSM_obj.data #Same dimension as data_arr
            
    def csv_record(self,para,sc,comb_data):
        #Method for recording data as a CSV file.
        start_dtstr = self.start_dt.strftime('%Y-%m-%dT%H-%M-%SZ')
        stop_dtstr = self.stop_dt.strftime('%Y-%m-%dT%H-%M-%SZ')
        csv_filename = '_'.join([para,sc,start_dtstr,stop_dtstr]) + '.csv'
        with open(csv_filename,'w', newline = '') as csvfile:
            csvwriter = csv.writer(csvfile)
                        
            for i in np.arange(len(comb_data[:,0])):
                csvwriter.writerow(comb_data[i,:])
                        
        csvfile.close()

"""
Example Script (Comment out if import module elsewhere)
"""
"""
start_dt = dt.datetime(2002,3,18,14,15,0)
stop_dt = dt.datetime(2002,3,18,15,15,0)
sclist = ['C1']
#sclist = ['C1','C2','C3','C4'] #for Fig. S2
paralist = ['ion_fluxU','e_fluxU','Bvec','ion_vel','SC_Pot'] #always put SCPot last, Fig 3
#paralist = ['H+_fluxU','He+_fluxU','O+_fluxU'] #Fig 4A-4C
#paralist = ['H+_fluxPA','He+_fluxPA','O+_fluxPA','e_fluxPA','e_fluxU180','SC_Pot'] #Fig 4E-4I
#paralist = ['e_fluxU','SC_Pot'] #Fig S2, all four SCs and sc_comp=True


instrlist = ['FGM', 'CIS', 'PEA', 'EFW']
saveformat = 'png'

testcall = CPanel(start_dt, stop_dt, sclist, paralist, instrlist, saveformat,
                  new_dl=True, veccomb=True, sc_comp=False, autoplot=True,
                  MP_HL=True, Walen_HL=True, label_start='B', xax_view=1, mode='light', 
                  data_record=True)
"""