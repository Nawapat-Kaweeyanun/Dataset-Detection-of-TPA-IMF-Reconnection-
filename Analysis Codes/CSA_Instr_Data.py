"""
Edition Date: 2024-August-08
@author: Nawapat Kaweeyanun
"""

"""
Objective: Create a class object for handling data extracted from one Cluster Science Archive CDF data file

Prerequisite: cdflib Python module

"""
import os
import shutil
import cdflib
import csv
import numpy as np
import datetime as dt
import CSA_Extract as CEx
import CSA_Organiser as COrg

"""
Class for Cluster Data

Instrument List (Cluster):
    FGM (fluxgate magnetometer)
    EFW (electric field double probe antenna)
    CIS (ion data)
    PEACE (electron data)
    RAPID (energetic ion and electron spectrometer)
    EDI (electron drift instrument)
    
Note: Not all parameters are available for certain instruments
"""

class Cdata(object):
    
    def __init__(self, cdffilepath):
        
        #determine data ID of CDF file
        cdffile = cdflib.CDF(cdffilepath)
        self.DataID = cdffilepath.split('/')[-1].split('__')[0]
        
        #determine which instrument, then call relevant extracting function
        if 'FGM' in self.DataID:
            self.FGM_Extract(cdffile)
        elif 'EFW' in self.DataID:
            self.EFW_Extract(cdffile)
        elif 'CIS' in self.DataID:
            self.CIS_Extract(cdffile)
        elif 'RAP' in self.DataID:
            self.RAP_Extract(cdffile)
        elif 'EDI' in self.DataID:
            self.EDI_Extract(cdffile)
        elif 'PEA' in self.DataID:
            self.PEA_Extract(cdffile)
        else:
            print('Extraction method unavailable for this instrument')
            return
        
    def FGM_Extract(self,cdffile):
        Epochdates = cdffile.varget('time_tags__' + self.DataID)
        self.datetime = np.array(cdflib.cdfepoch.to_datetime(Epochdates))
        self.hint = cdffile.varget('half_interval__' + self.DataID)
        self.Bvec = cdffile.varget('B_vec_xyz_gse__' + self.DataID)  
        self.Bmag = cdffile.varget('B_mag__' + self.DataID)
        self.SC_gse = cdffile.varget('sc_pos_xyz_gse__' + self.DataID)
        self.FGMrange = cdffile.varget('range__' + self.DataID)
        self.FGMtele = cdffile.varget('tm__' + self.DataID)
        self.periodname = self.periodname_calc()     
    
    def EFW_Extract(self,cdffile):
        Epochdates = cdffile.varget('time_tags__' + self.DataID)
        self.datetime = np.array(cdflib.cdfepoch.to_datetime(Epochdates))
        self.periodname = self.periodname_calc() 

        if 'L3_E' in self.DataID: #electric field        
            self.Evec = cdffile.varget('E_Vec_xyz_GSE__' + self.DataID)        
            self.deltaEz = cdffile.varget('delta_Ez_GSE__' + self.DataID)       
            self.Ebitmask = cdffile.varget('E_bitmask__' + self.DataID)        
            self.quality = cdffile.varget('E_quality__' + self.DataID)
            self.periodname = self.periodname_calc()  
        elif 'L3_P' in self.DataID: #spacecraft potential
            self.SC_Pot = cdffile.varget('Spacecraft_potential__' + self.DataID)
            self.probes = cdffile.varget('P_probes__' + self.DataID)
            self.aspoc = cdffile.varget('P_probes__' + self.DataID)
            self.bitmask = cdffile.varget('P_bitmask__' + self.DataID)
            self.quality = cdffile.varget('P_quality__' + self.DataID)
        
    def CIS_Extract(self,cdffile):    
        #check if the CIS ID is moments (bulk ions or species-specific)
        #or pitch angle particle differential flux distribution
        
        if 'MOMENTS' in self.DataID and 'HIA' in self.DataID: #ion moments
            Epochdates = cdffile.varget('time_tags__' + self.DataID)
            self.datetime = np.array(cdflib.cdfepoch.to_datetime(Epochdates))            
            self.hint = cdffile.varget('delta_time__' + self.DataID)
            self.sensitivity = cdffile.varget('sensitivity__' + self.DataID)
            self.CISmode = cdffile.varget('cis_mode__' + self.DataID)
            self.density = cdffile.varget('density__' + self.DataID)
            self.vel_isr2 = cdffile.varget('velocity_isr2__' + self.DataID)
            self.vel_gse = cdffile.varget('velocity_gse__' + self.DataID)
            self.temp = cdffile.varget('temperature__' + self.DataID)
            self.temp_Bll = cdffile.varget('temp_par__' + self.DataID)
            self.temp_Bperp = cdffile.varget('temp_perp__' + self.DataID)
            self.pres = cdffile.varget('pressure__' + self.DataID)
            self.prestens =  cdffile.varget('pressure_tensor__' + self.DataID)
            self.periodname = self.periodname_calc()  
                        
        elif (('MOMENTS' in self.DataID and 'H1' in self.DataID) or 
            ('MOMENTS' in self.DataID and 'He1' in self.DataID) or
            ('MOMENTS' in self.DataID and 'O1' in self.DataID)): #H+,He+,O+ moments
            Epochdates = cdffile.varget('time_tags__' + self.DataID)
            self.datetime = np.array(cdflib.cdfepoch.to_datetime(Epochdates))
            self.hint = cdffile.varget('duration__' + self.DataID)
            self.density = cdffile.varget('density__' + self.DataID)
            self.vel = cdffile.varget('velocity__' + self.DataID)
            self.pres = cdffile.varget('pressure__' + self.DataID)
            self.temp = cdffile.varget('T__' + self.DataID)
            self.temp_Bll = cdffile.varget('T_par__' + self.DataID)
            self.temp_Bperp = cdffile.varget('T_perp__' + self.DataID)
            self.periodname = self.periodname_calc()
                        
        elif 'PAD' in self.DataID: #differential flux distributions
            Epochdates = cdffile.varget('time_tags__' + self.DataID)
            self.datetime = np.array(cdflib.cdfepoch.to_datetime(Epochdates)) 
            self.hint = cdffile.varget('duration__' + self.DataID)
            self.CISmode = cdffile.varget('cis_mode__' + self.DataID)
            self.CIStele = cdffile.varget('tm_product__' + self.DataID)
            self.diffflux = cdffile.varget('Differential_Particle_Flux__' + self.DataID)
            self.pa = cdffile.varget('pitch_angle__' + self.DataID)
            self.energy = cdffile.varget('energy_table__' + self.DataID)
            self.energy_delplus = cdffile.varget('delta_plus_energy_table__' + self.DataID)
            self.energy_delminus = cdffile.varget('delta_minus_energy_table__' + self.DataID)
            self.periodname = self.periodname_calc()
        
        elif 'PF' in self.DataID and 'PAD' not in self.DataID: #3D pitch angle flux distributions
            Epochdates = cdffile.varget('time_tags__' + self.DataID)
            self.datetime = np.array(cdflib.cdfepoch.to_datetime(Epochdates)) 
            self.hint = cdffile.varget('duration__' + self.DataID)
            self.CISmode = cdffile.varget('cis_mode__' + self.DataID)
            self.CIStele = cdffile.varget('tm_product__' + self.DataID)
            self.diffflux = cdffile.varget('3d_ions__' + self.DataID)
            self.theta = cdffile.varget('theta__' + self.DataID)
            self.phi = cdffile.varget('phi__' + self.DataID)
            self.energy = cdffile.varget('energy_table__' + self.DataID)
            self.energy_delplus = cdffile.varget('delta_plus_energy_table__' + self.DataID)
            self.energy_delminus = cdffile.varget('delta_minus_energy_table__' + self.DataID)
            self.periodname = self.periodname_calc()
            
            
    def PEA_Extract(self,cdffile):
        #only cover pitch angle distributions, not 3D distributions
        if 'PITCH' in self.DataID:
            Epochdates = cdffile.varget('time_tags__' + self.DataID)
            self.datetime = np.array(cdflib.cdfepoch.to_datetime(Epochdates)) 
            self.diffflux = cdffile.varget('Data__' + self.DataID)
            self.diffflux_bg = cdffile.varget('BackgroundLevel__' + self.DataID)
            self.energy = cdffile.varget('Sweep_Energy__' + self.DataID) #per bin, not one list
            self.pa = cdffile.varget('Sweep_PitchAngle__' + self.DataID)
            self.quality = cdffile.varget('Status_Quality__' + self.DataID)
            self.periodname = self.periodname_calc()
        
        
    def RAP_Extract(self,cdffile): 
        #only cover omni-directional data
        Epochdates = cdffile.varget('Time_tags__' + self.DataID)
        self.datetime = np.array(cdflib.cdfepoch.to_datetime(Epochdates))
        self.hint = cdffile.varget('Time_half_interval__' + self.DataID)
        self.energy = cdffile.varget('Dimension_E__' + self.DataID)
        self.quality = cdffile.varget('Quality__' + self.DataID)
        self.periodname = self.periodname_calc()
        
        if 'ESPCT6' in self.DataID:
            self.diffflux = cdffile.varget('Electron_Dif_flux__' + self.DataID)
            self.diffflux_SD = cdffile.varget('Electron_Dif_flux_SD__' + self.DataID)
        elif 'He' in self.DataID:
            self.diffflux = cdffile.varget('Proton_Dif_flux__' + self.DataID)
            self.diffflux_SD = cdffile.varget('Proton_Dif_flux_SD__' + self.DataID)
        elif 'H' in self.DataID:
            self.diffflux = cdffile.varget('Helium_Dif_flux__' + self.DataID)
            self.diffflux_SD = cdffile.varget('Helium_Dif_flux_SD__' + self.DataID)

    def EDI_Extract(self,cdffile):
        Epochdates = cdffile.varget('time_tags__' + self.DataID)
        self.datetime = np.array(cdflib.cdfepoch.to_datetime(Epochdates)) 
        self.hint = cdffile.varget('half_interval__' + self.DataID)
        self.status = cdffile.varget('status__' + self.DataID)
        self.vel = cdffile.varget('V_ed_xyz_gse__' + self.DataID)
        self.Evec = cdffile.varget('E_xyz_gse__' + self.DataID)
        self.periodname = self.periodname_calc()
                
        if 'MP' in self.DataID: #Reduced Chi-sq only for best res
            self.Chi = cdffile.varget('Reduced_chi_sq__' + self.DataID)
     
    def periodname_calc(self):
        #set period name of for extracted data array variable
        if type(self.datetime[0]) is dt.datetime and type(self.datetime[-1]) is dt.datetime:
            start_datetime = self.datetime[0] - dt.timedelta(microseconds=self.datetime[0].microsecond) #remove microseconds
            end_datetime = self.datetime[-1] - dt.timedelta(microseconds=self.datetime[-1].microsecond)
            startname = start_datetime.strftime('%Y-%m-%dT%H-%M-%S')
            endname = end_datetime.strftime('%Y-%m-%dT%H-%M-%S')
            periodname = startname + "___" + endname
        else: #certain periods have no data and present empty arrays
            periodname = 'DTimeCorrupted'
        return periodname
    
    def csvwrite(self,savedir,sc):
        #Write extracted CDF data to a CSV file
        #Call as separate function after class initialisation (extraction)
        #Below is FGM only, but can adjust to different instruments
        os.makedirs(savedir,exist_ok=True)
        
        #generate csv filename
        datafilename = savedir + 'ClusterFGM_{}_'.format(sc) + self.periodname + ".csv"
        #print(datafilename)
        
        #define field names (column headers)
        fields = ['Datetime','Half_Interval','Bvec','Bmag','SC_gse','FGMRange','FGMTelemetry']
        
        #create csv.writer object
        with open(datafilename, 'w', newline = '') as csvfile: #newline remove blank rows
            csvwriter = csv.writer(csvfile)
             
            #write the fields first
            csvwriter.writerow(fields)
        
            #then write the columns
            for n in np.arange(len(self.datetime)):
                csvwriter.writerow([self.datetime[n],self.hint[n],self.Bvec[n],self.Bmag[n],
                                    self.SC_gse[n],self.FGMrange[n],self.FGMtele[n]])


"""
Script: Write CSV File for all data files selected (un-comment to use)
"""
"""
#Set up range of datetimes to write data
sc_choice = ['C1','C2','C3','C4']
PathStart = 'Datafiles_ClusterFGM'
SaveTarget = 'Datafiles_ClusterFGM_CSV'
yr_choice = ['2002']
mo_choice = [3]
day_choice = [18]
hr_range = [14,16] #[start,end] [inclusive,exclusive]

#File error list
problem_files = []

#If no spacecraft preference, trace all available spacecrafts.
if len(sc_choice) == 0:
    sclist = os.listdir(PathStart)
    sclist.sort()
else:
    sclist = sc_choice
    
#Loop over each spacecraft selected.    
for sc in sclist:
    
    #Obtain correct data ID.
    if sc == 'C1':
        DataID = 'C1_CP_FGM_SPIN'
    elif sc == 'C2':
        DataID = 'C2_CP_FGM_SPIN'
    elif sc == 'C3':
        DataID = 'C3_CP_FGM_SPIN'
    elif sc == 'C4':
        DataID = 'C4_CP_FGM_SPIN'
    
    #If no year preference, write all available years.
    if len(yr_choice) == 0:
        yr_list = os.listdir(PathStart + '{}/'.format(DataID))
        yr_list.sort()
        if '.DS_Store' in yr_list: #Remove .DS_Store file for Mac
            yr_list.remove('.DS_Store')
    else:
        yr_list = yr_choice 
    
    #Loop over each year.
    for yr in yr_list:
        
        #If no month preference, write all available months.
        if len(mo_choice) == 0:
            mo_list = os.listdir(PathStart + '{}/{}/'.format(DataID,yr))
            mo_list.sort()
            if '.DS_Store' in mo_list:
                mo_list.remove('.DS_Store')
        else:
            mo_list = mo_choice
        
        #Loop over each month.
        for mo in mo_list:
            #Obtain file list and sort alphabetically.
            FileFolderPath = PathStart + '/{}/{}/{}/'.format(DataID,yr,mo)
            SaveFolder = SaveTarget + '/{}/{}/{}/'.format(DataID,yr,mo)
            Filelist = os.listdir(FileFolderPath)
            Filelist.sort()
            
            #Apply day filter (if exist) and update file list
            if len(day_choice) != 0:
                Filelist2 = []
                for day in day_choice:
                    if day < 10 and mo < 10:
                        datestr = '{}0{}0{}'.format(yr,mo,day)
                    elif day < 10 and mo >= 10:
                        datestr = '{}{}0{}'.format(yr,mo,day)
                    elif day >= 10 and mo < 10:
                        datestr = '{}0{}{}'.format(yr,mo,day)
                    elif day >= 10 and mo >= 10:
                        datestr = '{}{}{}'.format(yr,mo,day)
                           
                    for file in Filelist: #Filter file with right date
                        if datestr in file.split('__')[1].split('_')[0]: #only look at start date
                            Filelist2.append(file)
                    
                Filelist2.sort()
                Filelist = Filelist2 #replace old file list with new one
            
            #Apply hour filter (if exist) and update file list    
            if len(hr_range) != 0:
                Filelist2 = []
                
                #Obtain filter hour strings.
                if hr_range[0] < 10:
                    hri_str = '0{}0000'.format(hr_range[0])
                else:
                    hri_str = '{}0000'.format(hr_range[0])
                if hr_range[1] < 10:
                    hrf_str = '0{}0000'.format(hr_range[1])
                elif hr_range[1] == 24:
                    hrf_str = 'None'
                else:
                    hrf_str = '{}0000'.format(hr_range[1])
                
                hri_indlist = []
                hrf_indlist = []
                
                #Find files that match filter hour strings.
                for m in np.arange(len(Filelist)):
                    if hri_str in Filelist[m].split('__')[1].split('_')[1]:
                        hri_indlist.append(m)
                        
                    if hrf_str == 'None':
                        hrf_indlist.append('None')
                    elif hrf_str in Filelist[m].split('__')[1].split('_')[1]:
                        hrf_indlist.append(m)
                        
                #If there is no file beyond specified time interval
                if len(hrf_indlist) == 0:
                    hrf_indlist.append('None')
                
                
                for n in np.arange(len(hri_indlist)):
                    
                    if hrf_indlist[n] == 'None':
                        Filelist2.append(Filelist[hri_indlist[n]:])
                    else:
                        #get files between start/end hours of each day
                        Filelist2.append(Filelist[hri_indlist[n]:hrf_indlist[n]])
                    
                Filelist2 = Filelist2[0] #remove extra array bracket
                Filelist2.sort()
                Filelist = Filelist2
            
            #Loop over all files to write
            for file_no in np.arange(len(Filelist)):

                Filename = Filelist[file_no]
                Filepath = FileFolderPath + '/' + Filename
                FileDataID = Filename.split('__')[0]
                FileStartDate = Filename.split('__')[1].split('_')[0]
                
                #Use try-except to skip error when there is no data
                try:
                    Para = Cdata(Filepath)
                    Para.csvwrite(SaveFolder,sc)
                    print('CSV file saved at {}-{}-{}, file no. {}'.format(sc,yr,mo,file_no+1))
                except Exception as e:
                    print(e)
                    print('Error prevented csv file writing at {}-{}-{}, file no.{}'.format(sc,yr,mo,file_no+1))
                    problem_files.append(Filename)
                    continue

"""

