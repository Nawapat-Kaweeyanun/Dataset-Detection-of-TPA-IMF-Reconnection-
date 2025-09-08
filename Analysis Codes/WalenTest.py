"""
Edition Date: 2025-September-01
@author: Nawapat Kaweeyanun
"""

"""
Objective: Perform Walen test for specified date time interval

Note: Script for application is provided at the end.
"""

import os
import shutil
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import CSA_Extract as CEx
import CSA_Organiser as COrg
import CSA_Instr_Data as CID
from scipy import optimize
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from matplotlib.ticker import AutoMinorLocator
import csv


class Walen(object):

    def __init__(self,SC,test_dti,test_dtf,ref_dti=None,ref_dtf=None,
                 sub_only=False,title=True,mode='light',saveformat='png',data_record=False):
        self.SC = SC #only one spacecraft permitted
        self.test_dti = test_dti
        self.test_dtf = test_dtf
        self.ref_dti = ref_dti #Note: reference time period needed for CSA download
        self.ref_dtf = ref_dtf

        #Set Cluster instrument list required
        instrlist = ['FGM','CIS']

        for instr in instrlist:
                    
			#Download CIS data (see method for details) & allocate to FGM or CIS
            TestData = self.CSADownload(self.test_dti,self.test_dtf,self.SC,instr)
            if instr == 'FGM':
                self.FGMTestData = TestData
            elif instr == 'CIS':
                self.CISTestData = TestData
                    
        #If density substitute is applied i.e., reference datetime exists, also download relevant information            
        if self.ref_dti is not None:
                            
            for instr in instrlist:
                RefData = self.CSADownload(self.ref_dti,self.ref_dtf,self.SC,instr)
                if instr == 'FGM':
                    self.FGMRefData = RefData
                elif instr == 'CIS':
                    self.CISRefData = RefData
                        
                        
        #Extract magnetic field data in test interval from FGM, then convert from GSE to GSM.
        self.FGMTestdatetime = self.FGMTestData.datetime
        self.FGMTestdatetime = self.FGMTestData.datetime
        BvecTestGSM = self.GSM_from_GSE(self.FGMTestdatetime,self.FGMTestData.Bvec)
        self.BxTest = BvecTestGSM[:,0]
        self.ByTest = BvecTestGSM[:,1]
        self.BzTest = BvecTestGSM[:,2]
        
        #Extract ion moments data from CIS in test interval.        
        self.CISTestdatetime = self.CISTestData.datetime
        self.TllTest = self.CISTestData.temp_Bll
        self.TperpTest = self.CISTestData.temp_Bperp
        self.nTest = self.CISTestData.density
        
        #Convert V Test data from CIS to GSM
        VvecTestGSM = self.GSM_from_GSE(self.CISTestdatetime,self.CISTestData.vel_gse)
        self.vxTest = VvecTestGSM[:,0] #does not work if there is only one data point
        self.vyTest = VvecTestGSM[:,1]
        self.vzTest = VvecTestGSM[:,2]
        
        if self.ref_dti is not None:
            
            #Convert magnetic field and velocity data in the reference time range from GSE to GSM.
            BvecRefGSM = self.GSM_from_GSE(self.FGMRefData.datetime,self.FGMRefData.Bvec)
            VvecRefGSM = self.GSM_from_GSE(self.CISRefData.datetime,self.CISRefData.vel_gse)
    
            #Extract magnetic field data at referent point - take first datetime of reference period only
            self.FGMRefdatetime = self.FGMRefData.datetime[0]
            self.BxRef = BvecRefGSM[0,0]
            self.ByRef = BvecRefGSM[0,1]
            self.BzRef = BvecRefGSM[0,2]
    
           #Extract temperature, pressure, and velocity data at referent point - first datetime only.
            self.CISRefdatetime = self.CISRefData.datetime[0]
            self.TllRef = self.CISRefData.temp_Bll[0]
            self.TperpRef = self.CISRefData.temp_Bperp[0]
            self.nRef = self.CISRefData.density[0]
            self.vxRef = VvecRefGSM[0,0]
            self.vyRef = VvecRefGSM[0,1]
            self.vzRef = VvecRefGSM[0,2]
        

        #Check that FGM Test and CIS Test datetime arrays have the same lengths
        #Implement averaging procedure if not match
        if len(self.FGMTestdatetime) != len(self.CISTestdatetime):
            #Can be gap in B-field. Rob interpolated.
            print("Datetime lengths between FGM and CIS not matched.") #reinstate raise once sort out DS averaging
            
            #Remove microseconds
            FGM_dt_no_ms = [dtime-dt.timedelta(microseconds=dtime.microsecond) for dtime in self.FGMTestdatetime]
            CIS_dt_no_ms = [dtime-dt.timedelta(microseconds=dtime.microsecond) for dtime in self.CISTestdatetime]
            
            #Find datetimes common between the two lists
            dt_matches = list(set(FGM_dt_no_ms).intersection(CIS_dt_no_ms))
            dt_matches = np.array(dt_matches)
            dt_matches.sort() #making sure datetime list is in order, otherwise mess up future appending
            
            #Set up arrays to record parameters at common datetimes
            newdt = []
            newBx = []
            newBy = []
            newBz = []
            newTll = []
            newTperp = []
            newn = []
            newvx = []
            newvy = []
            newvz = []
            
            #Loop over each common datetime
            for dtime in dt_matches:
                dtime_F = np.array([dtime for FGMdt in self.FGMTestdatetime])
                dtime_C = np.array([dtime for CISdt in self.CISTestdatetime])
                FGM_diffarr = self.FGMTestdatetime - dtime_F #same second
                CIS_diffarr = self.CISTestdatetime - dtime_C
                
                FGMbool = []
                for item in FGM_diffarr:
                    if item.seconds == 0:
                        FGMbool.append(True)
                    else:
                        FGMbool.append(False)
                FGMbool = np.array(FGMbool)
                
                CISbool = []
                for item in CIS_diffarr:
                    if item.seconds==0:
                        CISbool.append(True)
                    else:
                        CISbool.append(False)
                CISbool = np.array(CISbool)
                
                #take first element to remove array + useful for multiple matches. Use FGM for dt. Doesn't matter which is longer.
                newdt.append(self.FGMTestdatetime[FGMbool][0]) 
                newBx.append(self.BxTest[FGMbool][0])
                newBy.append(self.ByTest[FGMbool][0])
                newBz.append(self.BzTest[FGMbool][0])
                newTll.append(self.TllTest[CISbool][0])
                newTperp.append(self.TperpTest[CISbool][0])
                newn.append(self.nTest[CISbool][0])
                newvx.append(self.vxTest[CISbool][0])
                newvy.append(self.vyTest[CISbool][0])
                newvz.append(self.vzTest[CISbool][0])
                
            #Convert list to array    
            self.FGMTestdatetime = np.array(newdt)
            self.CISTestdatetime = np.array(newdt)
            self.BxTest = np.array(newBx)
            self.ByTest = np.array(newBy)            
            self.BzTest = np.array(newBz)
            self.TllTest = np.array(newTll)
            self.TperpTest = np.array(newTperp)
            self.nTest = np.array(newn)
            self.vxTest = np.array(newvx)
            self.vyTest = np.array(newvy)            
            self.vzTest = np.array(newvz)
            print('Datetime arrays length sorted')

        
        #Calculate magnetopause normal vector (just in case)
        self.MPn = self.MVAB(self.BxTest,self.ByTest,self.BzTest)

        #Find time in seconds (with decimals) between plotting datetimes from the first element of datetime list
        delta_arr = self.FGMTestdatetime - self.FGMTestdatetime[0]
        self.t_arr = np.array([delta.seconds+delta.microseconds/10**6 for delta in delta_arr])
        
        #Calculate ion velocities in de Hoffman-Teller frame (see method for details)
        self.vHTTest, self.aHT, self.vHT0 = self.vHTCalc(self.t_arr,self.BxTest,
                                                         self.ByTest,self.BzTest,
                                                         self.vxTest,self.vyTest,
                                                         self.vzTest)
    
        #Calculate convective and HT electric field (in V/m) (see method for details)
        self.ECxTest,self.ECyTest,self.ECzTest = self.ECalc(self.BxTest*10**-9,self.ByTest*10**-9,self.BzTest*10**-9,
                                                            self.vxTest*10**3,self.vyTest*10**3,self.vzTest*10**3)
        self.EHTxTest,self.EHTyTest,self.EHTzTest = self.ECalc(self.BxTest*10**-9,self.ByTest*10**-9,self.BzTest*10**-9,
                                                               self.vHTTest[:,0]*10**3,self.vHTTest[:,1]*10**3,
                                                               self.vHTTest[:,2]*10**3)
        self.ECxTest = self.ECxTest*10**3 #convert from V/m to mV/m
        self.ECyTest = self.ECyTest*10**3 
        self.ECzTest = self.ECzTest*10**3 
        self.EHTxTest = self.EHTxTest*10**3
        self.EHTyTest = self.EHTyTest*10**3 
        self.EHTzTest = self.EHTzTest*10**3 
        
        #Find ion velocity in HT frame.
        self.vdiffxTest = self.vxTest - self.vHTTest[:,0]
        self.vdiffyTest = self.vyTest - self.vHTTest[:,1]
        self.vdiffzTest = self.vzTest - self.vHTTest[:,2]

        #Calculate Alfven velocity in HT frame (see method for details)
        if self.ref_dti is None:
            self.vAxTest,self.vAyTest,self.vAzTest = self.vACalc(self.BxTest,self.ByTest,
                                                                 self.BzTest,self.nTest,
                                                                 self.TllTest,self.TperpTest)
        else:
            self.vAxTest,self.vAyTest,self.vAzTest,self.vAxSub,self.vAySub,self.vAzSub = \
                self.vACalc(self.BxTest,self.ByTest,
                            self.BzTest,self.nTest,
                            self.TllTest,self.TperpTest,
                            self.BxRef,self.ByRef,
                            self.BzRef,self.nRef,
                            self.TllRef,self.TperpRef)


        #Define figure properties for light and dark mode
        if mode == 'light':
            ctick = 'k' #Black ticks/labels
            tchoice = False #Plot with white background
        elif mode == 'dark':
            ctick = 'w' #White ticks/labels
            tchoice = True #Plot with transparent background
        else: #raise exception for incorrect input
            raise Exception('Incorrect Colour Mode')

        #Call plotting function
        #Consider three cases (1) 2-plots, no reference-point density substitution,
        #(2) 2-plots with substituion, and (3) 3-plots with substitution
        if self.ref_dti is None:
            fig,ax = plt.subplots(1,2,figsize=(12,6))
            fig,a0,b0,corr0 = self.WalenPlot(fig,ax[0],mode,self.EHTxTest,self.EHTyTest,self.EHTzTest,
                           self.ECxTest,self.ECyTest,self.ECzTest)
            fig,a1,b1,corr1 = self.WalenPlot(fig,ax[1],mode,self.vAxTest,self.vAyTest,self.vAzTest,
                           self.vdiffxTest,self.vdiffyTest,self.vdiffzTest)
            ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
            ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))
            ax[0].tick_params(axis='both',which='major',direction='inout',width=1,length=15)
            ax[0].tick_params(axis='both',which='minor',direction='inout',width=1,length=5)
            ax[0].tick_params(axis='both', which='major', labelsize=14)
            ax[0].set_xlabel('$\mathrm{E_{HT,\gamma}}$ (V/m)',fontsize=16,fontfamily='sans-serif')
            ax[0].set_ylabel('$\mathrm{E_{C,\gamma}}$ (V/m)',fontsize=16,fontfamily='sans-serif')
            ax[0].axis('equal')
            ax[0].title.set_text('(a) Corr. {}, y-int. {}'.format(np.round(corr0,2),np.round(b0,2)))
            ax[0].title.set_fontsize(16)
            ax[0].title.set_fontfamily('sans-serif')
            ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))
            ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))
            ax[1].tick_params(axis='both',which='major',direction='inout',width=1,length=15)
            ax[1].tick_params(axis='both',which='minor',direction='inout',width=1,length=5)
            ax[1].tick_params(axis='both', which='major', labelsize=14)
            ax[1].set_xlabel('$\mathrm{v_{ATest,\gamma}}$ (km/s)',fontsize=16,fontfamily='sans-serif')
            ax[1].set_ylabel('$\mathrm{(v-v_{HT})_\gamma}$ (km/s)',fontsize=16,fontfamily='sans-serif')  
            ax[1].axis('equal')
            ax[1].title.set_text('(b) Corr. {}, y-int. {}'.format(np.round(corr1,2),np.round(b1,2)))
            ax[1].title.set_fontsize(16)
            ax[1].title.set_fontfamily('sans-serif')
            self.fit_grad = [a0,a1] #calculate gradients, y-intercepts, and correlation coefficients 
            self.fit_yint = [a0,a1]
            self.corr_coeff = [corr0,corr1]
            
            #Set colors
            ax[0].spines['bottom'].set_color(ctick) #Set axes color
            ax[0].spines['top'].set_color(ctick)
            ax[0].spines['left'].set_color(ctick)
            ax[0].spines['right'].set_color(ctick)
            ax[0].tick_params(axis='both',which='both',colors=ctick) #Set tick color
            ax[0].xaxis.label.set_color(ctick) #Set x-axis ticklabel color
            ax[0].yaxis.label.set_color(ctick) #Set y-axis ticklabel color
            ax[0].title.set_color(ctick)
            ax[1].spines['bottom'].set_color(ctick) #Set axes color
            ax[1].spines['top'].set_color(ctick)
            ax[1].spines['left'].set_color(ctick)
            ax[1].spines['right'].set_color(ctick)
            ax[1].tick_params(axis='both',which='both',colors=ctick) #Set tick color
            ax[1].xaxis.label.set_color(ctick) #Set x-axis ticklabel color
            ax[1].yaxis.label.set_color(ctick) #Set y-axis ticklabel color
            ax[1].title.set_color(ctick)
            
            #If required, save plot data to CSV
            if data_record == True:
                self.csv_recorder2()
                print('2-Panel Walen data recorded in csv')
        
        elif self.ref_dti is not None and sub_only == True:
            fig,ax = plt.subplots(1,2,figsize=(12,6))
            fig,a0,b0,corr0 = self.WalenPlot(fig,ax[0],mode,self.EHTxTest,self.EHTyTest,self.EHTzTest,
                           self.ECxTest,self.ECyTest,self.ECzTest)
            fig,a1,b1,corr1 = self.WalenPlot(fig,ax[1],mode,self.vAxSub,self.vAySub,self.vAzSub,
                                       self.vdiffxTest,self.vdiffyTest,self.vdiffzTest)
            ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
            ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))
            ax[0].tick_params(axis='both',which='major',direction='inout',width=1,length=15)
            ax[0].tick_params(axis='both',which='minor',direction='inout',width=1,length=5)
            ax[0].tick_params(axis='both', which='major', labelsize=14)
            ax[0].set_xlabel('$\mathrm{E_{HT,\gamma}}$ (V/m)',fontsize=16,fontfamily='sans-serif')
            ax[0].set_ylabel('$\mathrm{E_{C,\gamma}}$ (V/m)',fontsize=16,fontfamily='sans-serif')
            ax[0].axis('equal')
            ax[0].title.set_text('(a) Corr. {}, y-int. {}'.format(np.round(corr0,2),np.round(b0,2)))
            ax[0].title.set_fontsize(16)
            ax[0].title.set_fontfamily('sans-serif')
            ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))
            ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))
            ax[1].tick_params(axis='both',which='major',direction='inout',width=1,length=15)
            ax[1].tick_params(axis='both',which='minor',direction='inout',width=1,length=5)
            ax[1].tick_params(axis='both', which='major', labelsize=14)
            ax[1].set_xlabel('$\mathrm{v_{ASub,\gamma}}$ (km/s)',fontsize=16,fontfamily='sans-serif')
            ax[1].set_ylabel('$\mathrm{(v-v_{HT})_\gamma}$ (km/s)',fontsize=16,fontfamily='sans-serif')
            ax[1].axis('equal')
            ax[1].title.set_text('(b) Corr. {}, y-int. {}'.format(np.round(corr1,2),np.round(b1,2)))
            ax[1].title.set_fontsize(16)
            ax[1].title.set_fontfamily('sans-serif')
            self.fit_grad = [a0,a1]
            self.fit_yint = [a0,a1]
            self.corr_coeff = [corr0,corr1]
            
            ax[0].spines['bottom'].set_color(ctick)
            ax[0].spines['top'].set_color(ctick)
            ax[0].spines['left'].set_color(ctick)
            ax[0].spines['right'].set_color(ctick)
            ax[0].tick_params(axis='both',which='both',colors=ctick) 
            ax[0].xaxis.label.set_color(ctick) 
            ax[0].yaxis.label.set_color(ctick) 
            ax[0].title.set_color(ctick)
            ax[1].spines['bottom'].set_color(ctick) 
            ax[1].spines['top'].set_color(ctick)
            ax[1].spines['left'].set_color(ctick)
            ax[1].spines['right'].set_color(ctick)
            ax[1].tick_params(axis='both',which='both',colors=ctick) 
            ax[1].xaxis.label.set_color(ctick) 
            ax[1].yaxis.label.set_color(ctick) 
            ax[1].title.set_color(ctick)
            
            if data_record == True:
                self.csv_recorder2()
                print('2-Panel Walen data recorded in csv')
            
        elif self.ref_dti is not None and sub_only == False:
            fig,ax = plt.subplots(1,3,figsize=(18,6))
            fig,a0,b0,corr0 = self.WalenPlot(fig,ax[0],mode,self.EHTxTest,self.EHTyTest,self.EHTzTest,
                           self.ECxTest,self.ECyTest,self.ECzTest)
            fig,a1,b1,corr1 = self.WalenPlot(fig,ax[1],mode,self.vAxTest,self.vAyTest,self.vAzTest,
                           self.vdiffxTest,self.vdiffyTest,self.vdiffzTest)
            fig,a2,b2,corr2 = self.WalenPlot(fig,ax[2],mode,self.vAxSub,self.vAySub,self.vAzSub,
                                       self.vdiffxTest,self.vdiffyTest,self.vdiffzTest)
            ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
            ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))
            ax[0].tick_params(axis='both',which='major',direction='inout',width=1,length=15)
            ax[0].tick_params(axis='both',which='minor',direction='inout',width=1,length=5)
            ax[0].tick_params(axis='both', which='major', labelsize=14)
            ax[0].set_xlabel('$\mathrm{E_{HT,\gamma}}$ (V/m)',fontsize=16,fontfamily='sans-serif')
            ax[0].set_ylabel('$\mathrm{E_{C,\gamma}}$ (V/m)',fontsize=16,fontfamily='sans-serif')
            ax[0].axis('equal')
            ax[0].title.set_text('(a) Corr. {}, y-int. {}'.format(np.round(corr0,2),np.round(b0,2)))
            ax[0].title.set_fontsize(16)
            ax[0].title.set_fontfamily('sans-serif')
            ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))
            ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))
            ax[1].tick_params(axis='both',which='major',direction='inout',width=1,length=15)
            ax[1].tick_params(axis='both',which='minor',direction='inout',width=1,length=5)
            ax[1].tick_params(axis='both', which='major', labelsize=14)
            ax[1].set_xlabel('$\mathrm{v_{ATest,\gamma}}$ (km/s)',fontsize=16,fontfamily='sans-serif')
            ax[1].set_ylabel('$\mathrm{(v-v_{HT})_\gamma}$ (km/s)',fontsize=16,fontfamily='sans-serif')  
            ax[1].axis('equal')
            ax[1].title.set_text('(b) Corr. {}, y-int. {}'.format(np.round(corr1,2),np.round(b1,2)))
            ax[1].title.set_fontsize(16)
            ax[1].title.set_fontfamily('sans-serif')
            ax[2].xaxis.set_minor_locator(AutoMinorLocator(5))
            ax[2].yaxis.set_minor_locator(AutoMinorLocator(5))
            ax[2].tick_params(axis='both',which='major',direction='inout',width=1,length=15)
            ax[2].tick_params(axis='both',which='minor',direction='inout',width=1,length=5)
            ax[2].tick_params(axis='both', which='major', labelsize=14)
            ax[2].set_xlabel('$\mathrm{v_{ASub,\gamma}}$ (km/s)',fontsize=16,fontfamily='sans-serif')
            ax[2].set_ylabel('$\mathrm{(v-v_{HT})_\gamma}$ (km/s)',fontsize=16,fontfamily='sans-serif')
            ax[2].axis('equal')
            ax[2].title.set_text('(c) Corr. {}, y-int. {}'.format(np.round(corr2,2),np.round(b2,2)))
            ax[2].title.set_fontsize(16)
            ax[2].title.set_fontfamily('sans-serif')
            self.fit_grad = [a0,a1,a2] 
            self.fit_yint = [b0,b1,b2]
            self.corr_coeff = [corr0,corr1,corr2]

            ax[0].spines['bottom'].set_color(ctick)
            ax[0].spines['top'].set_color(ctick)
            ax[0].spines['left'].set_color(ctick)
            ax[0].spines['right'].set_color(ctick)
            ax[0].tick_params(axis='both',which='both',colors=ctick)
            ax[0].xaxis.label.set_color(ctick)
            ax[0].yaxis.label.set_color(ctick)
            ax[0].title.set_color(ctick)
            ax[1].spines['bottom'].set_color(ctick)
            ax[1].spines['top'].set_color(ctick)
            ax[1].spines['left'].set_color(ctick)
            ax[1].spines['right'].set_color(ctick)
            ax[1].tick_params(axis='both',which='both',colors=ctick)
            ax[1].xaxis.label.set_color(ctick) 
            ax[1].yaxis.label.set_color(ctick) 
            ax[1].title.set_color(ctick)
            ax[2].spines['bottom'].set_color(ctick) 
            ax[2].spines['top'].set_color(ctick)
            ax[2].spines['left'].set_color(ctick)
            ax[2].spines['right'].set_color(ctick)
            ax[2].tick_params(axis='both',which='both',colors=ctick) 
            ax[2].xaxis.label.set_color(ctick) 
            ax[2].yaxis.label.set_color(ctick)
            ax[2].title.set_color(ctick)
            
            if data_record == True:
                self.csv_recorder3()
                print('3-Panel Walen data recorded in csv')
       
        #Set initial and final datetime strings for saving
        dti_str = self.test_dti.strftime('%Y-%m-%dT%H-%M-%SZ')
        dtf_str = self.test_dtf.strftime('%Y-%m-%dT%H-%M-%SZ')
        
        #Adjust space between subplots
        plt.subplots_adjust(left=0.08,right=0.97,bottom=0.14,top=0.92,wspace=0.3,
                            hspace=0.08)
        
        #Set plot title if required
        if self.ref_dti is None and title == True:
            plt.suptitle('Walen Test Interval:{}, {}-{}'.format(self.test_dti.strftime('%Y-%m-%d'),
                                                                self.test_dti.strftime('%H:%M:%S'),
                                                                self.test_dtf.strftime('%H:%M:%S')),
                         fontsize=18,color=ctick)
            
        elif self.ref_dti is not None and title == True:
            plt.suptitle('Walen Test Interval: {}, {}-{} \n Density Ref Time: {}'.format(self.test_dti.strftime('%Y-%m-%d'),
                                                                                       self.test_dti.strftime('%H:%M:%S'),
                                                                                       self.test_dtf.strftime('%H:%M:%S'),
                                                                                       self.ref_dti.strftime('%H:%M:%S')),
                         fontsize=18,color=ctick)
            
        #Save figure
        savename = 'WalenTest_{}_{}_{}.{}'.format(self.SC,dti_str,dtf_str,saveformat)
        plt.savefig(savename,dpi=100,format=saveformat,transparent=tchoice)

    def CSADownload(self,start_dt,stop_dt,SC,instr):
        #Method for downloading Walen-relevant data from Cluster Science Archive
        #See PanelPlotter.py for more details
        
        #Set up URL and datetimes
        myurl = 'https://csa.esac.esa.int/csa-sl-tap/data'
        start_dtstr = start_dt.strftime('%Y-%m-%dT%H:%M:%SZ')
        stop_dtstr = stop_dt.strftime('%Y-%m-%dT%H:%M:%SZ')
        
        #Remove previous run if exist
        if 'CSA_ToMove' in os.listdir(os.getcwd()): 
            shutil.rmtree('CSA_ToMove')
            os.makedirs('CSA_ToMove',exist_ok=True)

        #Obtain data ID
        if instr == 'FGM':
       		ID = '{}_CP_FGM_SPIN'.format(SC.upper())
       	elif instr == 'CIS' and 'C' in SC.upper(): #Cluster ion
       		ID = '{}_CP_CIS-HIA_ONBOARD_MOMENTS'.format(SC.upper())

        #Set parameter dictionary
       	query_specs = {'RETRIEVAL_TYPE': 'product',
                'DATASET_ID': ID,
                'START_DATE': start_dtstr,
                'END_DATE': stop_dtstr,
                'DELIVERY_FORMAT': 'CDF'}
          
        #Download data
        try:
        	CEx.download(myurl, query_specs, 'tap_download.tgz')
        	filepath = CEx.unpack_tgz('tap_download.tgz')
        except Exception:
            print(Exception)
            print('Error downloading {}. Check if data exist'.format(ID))

        #Move files to intermediate holding folder
        filedir = '/'.join(filepath.split('/')[0:-1])
        shutil.move(filedir,'CSA_ToMove' + '/' + filedir)
        
        #Set up final data folder
        targetfolder = 'ClusterFolder'
        if targetfolder in os.listdir(os.getcwd()):
            shutil.rmtree(targetfolder) #Remove previous run if exist.
        os.makedirs(targetfolder,exist_ok=True)     
        
        #Shift files to final folder
        COrg.fileshifter('CSA_ToMove',targetfolder)
        IDlist = os.listdir(targetfolder)
        
        #Remove CSA Download Folders to reduce clutter
        DownloadFolderList = [d for d in os.listdir(os.getcwd()) if 'CSA_Download' in d]
        for DLF in DownloadFolderList:
            try:
                shutil.rmtree(DLF) #Remove dir and files inside
            except: #Skip problematic/corrupted files.
                continue
        
        #Note: cannot append data across multiple Cluster data files (not an issue for short analysis periods)
        if len(IDlist) > 1:
            raise ('More than one file to extract. Not suited for this function.')

        #Find year and month (at period start) from saved data file path.
        yr = os.listdir('{}/{}'.format(targetfolder,IDlist[0]))[0]
        mo = os.listdir('{}/{}/{}'.format(targetfolder,ID,yr))[0]

        #Still can only work with one file
        filelist = os.listdir('{}/{}/{}/{}'.format(targetfolder,ID,yr,mo))
        if len(filelist) > 1:
            raise ('More than one file to extract. Not suited for this function.')
        
        #Construct final file path
        filepath = '{}/{}/{}/{}/{}'.format(targetfolder,ID,yr,mo,filelist[0])
        
        #Obtain Cluster data
        Data = CID.Cdata(filepath)
        
        return Data
    
    def GSM_from_GSE(self,dt_arr,data_arr):
        #Convert GSE data array to GSM
        #Take data array with datetime along rows.
        #Create SpacePy GSE 3D vector object and tag it with datetime, then convert it to GSM.

        GSE_obj = coord.Coords(data_arr,'GSE','car')
        GSE_obj.ticks = Ticktock(dt_arr,'UTC')
        GSM_obj = GSE_obj.convert('GSM','car')
        
        return GSM_obj.data #same dimension as data_arr


    def MVAB(self,Bx_arr,By_arr,Bz_arr):
        #Use minimum variance analysis to compute direction of MP normal
        
        #Construct magnetic variance matrix
        MB = np.zeros((3,3))
        
        M = len(Bx_arr)

        #Fill each element of variance matrix according to formula
        for i in np.arange(3): #Each component of GSE coordinates
            for j in np.arange(3):
            
                Bi_ave = 0
                Bj_ave = 0
                BiBj_ave = 0
                
                if i == 0:
                    Bi_arr = Bx_arr
                elif i == 1:
                    Bi_arr = By_arr
                elif i == 2:
                    Bi_arr = Bz_arr
                
                if j == 0:
                    Bj_arr = Bx_arr
                elif j == 1:
                    Bj_arr = By_arr
                elif j == 2:
                    Bj_arr = Bz_arr
                
                
                for m in np.arange(M):
                    
                    Bi = Bi_arr[m]
                    Bj = Bj_arr[m]
                    BiBj = Bi_arr[m]*Bj_arr[m]
                    
                    Bi_ave = Bi_ave + Bi
                    Bj_ave = Bj_ave + Bj
                    BiBj_ave = BiBj_ave + BiBj
                    
                
                Bi_ave = Bi_ave/M
                Bj_ave = Bj_ave/M
                BiBj_ave = BiBj_ave/M
                
                
                MB[i,j] = BiBj_ave - (Bi_ave*Bj_ave)
                
        
        #Find eigenvalues and eigenvectors
        eigenvalues, eigenvectors = np.linalg.eig(MB) #Does not always give eigenvalues/vectors in decreasing order
        
        #Find minimum of eigenvalue and its index
        ei_min = eigenvalues.min()
        ei_min_ind = np.where(eigenvalues == ei_min)[0][0]
        nvec = eigenvectors[ei_min_ind]
        
        return nvec
        
    
    def ECalc(self,Bx_arr,By_arr,Bz_arr,vx_arr,vy_arr,vz_arr):
        #Cross product to calculate electric field from magnetic field and velocity
        
        Ex_arr = -1*(vy_arr*Bz_arr - vz_arr*By_arr)
        Ey_arr = -1*(vz_arr*Bx_arr - vx_arr*Bz_arr)
        Ez_arr = -1*(vx_arr*By_arr - vy_arr*Bx_arr)
        
        return Ex_arr,Ey_arr,Ez_arr
    
    def vHTCalc(self,t_arr,Bx_arr,By_arr,Bz_arr,vx_arr,vy_arr,vz_arr):
        #Method for calculating velocity of de Hoffman-Teller frame
        
        #Set up Kronecker delta matrix
        kron_delta = np.array([[1,0,0],[0,1,0],[0,0,1]])
        
        #Loop over element then sum for K0, K1, and K2
        K0 = np.zeros((3,3)) #Zero 3x3 matrix
        K1 = np.zeros((3,3))
        K2 = np.zeros((3,3))
        Kv = np.zeros((3,1))
        Kvt = np.zeros((3,1))
        
        #Loop over all data points (use x-component since all three should be equal in length)
        M = len(Bx_arr)
        for i in np.arange(M):
            
            #Set up B matrices
            B_vector = np.array([Bx_arr[i],By_arr[i],Bz_arr[i]]) #Make 1D vector
            B_vert = B_vector[:,None] #Make 2D vertical column vector
            B_hor = B_vert.transpose() #Make 2D horizontal row vector
            Bsq_mat = np.matmul(B_vert,B_hor) #If inverted will get single element!
            Bsq_mag = Bx_arr[i]**2 + By_arr[i]**2 + Bz_arr[i]**2
            
            #Find this iteration of K0, K1, and K2
            K0_iter = Bsq_mag*kron_delta - Bsq_mat #K0_iter is K proj matrix itself 
            K1_iter = K0_iter*t_arr[i] 
            K2_iter = K0_iter*(t_arr[i])**2 
            
            #Sum
            K0 = K0 + K0_iter
            K1 = K1 + K1_iter
            K2 = K2 + K2_iter
            
            #Set up velocity vectors
            v_vector = np.array([vx_arr[i],vy_arr[i],vz_arr[i]]) #make 1D vector
            v_vert = v_vector[:,None] #make 2D column vector
            
            #Find this iteration of Kv, Kvt
            Kv_iter = np.matmul(K0_iter,v_vert) #get vertical column 
            Kvt_iter = Kv_iter*t_arr[i] 
            
            #Sum
            Kv = Kv + Kv_iter
            Kvt = Kvt + Kvt_iter
            
        #Average by M at the end of the summation
        K0 = K0/M
        K1 = K1/M
        K2 = K2/M
        Kv = Kv/M
        Kvt = Kvt/M
        
        
        #Compute acceleration and initial velocity of HT frame
        K0_inv = np.linalg.inv(K0) 
        K1_inv = np.linalg.inv(K1) 
        
        aHT_1st = np.linalg.inv(np.matmul(K0_inv,K1)-np.matmul(K1_inv,K2)) #First bracket of aHT expression
        aHT_2nd = np.matmul(K0_inv,Kv) - np.matmul(K1_inv,Kvt)
        aHT = np.matmul(aHT_1st,aHT_2nd) #Should get 2D vertical vector #Unit = km/s**2
    
        vHT0 = np.matmul(K0_inv,Kv - np.matmul(K1,aHT)) #should get 2D vertical vector #Unit = km/s
        
        #Convert back to length-3 1D array
        aHT = aHT[:,0]
        vHT0 = vHT0[:,0]
        
        #Obtain array of HT velocities (Mx3 array)
        vHT_arr = np.zeros((M,3))
        for m in np.arange(M):
            for j in np.arange(3):
                vHT_arr[m,j] = vHT0[j] + aHT[j]*t_arr[m]
        
        return vHT_arr,aHT,vHT0
    
    
    def vACalc(self,BxTest,ByTest,BzTest,nTest,TllTest,TperpTest,
               BxRef=None,ByRef=None,BzRef=None,nRef=None,TllRef=None,
               TperpRef=None):
        #Method to calcualte Alfven velocity in de Hoffman-Teller frame
        
        #Parameter set-up
        mu0 = 4*np.pi*10**-7
        kB = 1.38 * 10**-23
        m = 1.12*1.67*10**-27 #average ion mass, assumed 96% H+, 4% He2+
        PllTest = nTest*10**6*kB*TllTest/10**6 #in Pa, for 2 d.o.f. (para/perp)
        PperpTest = nTest*10**6*kB*TperpTest/10**6
        
        #Calculate mass density from number density
        rhoTest = m*nTest*10**6 #in kg/m^-3
        
        #Calculate pressure anisotropy factor
        BsqTest = (BxTest**2 + ByTest**2 + BzTest**2)*10**-18 #in T
        alphaTest = (PllTest - PperpTest)*mu0/BsqTest #dimensionless
        
        #Calculate Alfven velocity in km/s
        vAxTest = BxTest*10**-9*np.sqrt((1-alphaTest)/(mu0*rhoTest))
        vAyTest = ByTest*10**-9*np.sqrt((1-alphaTest)/(mu0*rhoTest))
        vAzTest = BzTest*10**-9*np.sqrt((1-alphaTest)/(mu0*rhoTest))
        vAxTest = vAxTest/1000
        vAyTest = vAyTest/1000
        vAzTest = vAzTest/1000
        
        #If substitution deployed, recalculate Alfven velocity
        if BxRef is not None: 
            PllRef = nRef*10**6*kB*TllRef/10**6
            PperpRef = nRef*10**6*kB*TperpRef/10**6
            
            rhoRef = m*nRef*10**6
            
            BsqRef = (BxRef**2 + ByRef**2 + BzRef**2)*10**-18 #in T
            alphaRef = (PllRef - PperpRef)*mu0/BsqRef #dimensionless
            
            rhoSub = rhoRef*(1-alphaRef)/(1-alphaTest)
            
            vAxSub = BxTest*10**-9*np.sqrt((1-alphaRef)/(mu0*rhoSub)) #use Btest, not sub
            vAySub = ByTest*10**-9*np.sqrt((1-alphaRef)/(mu0*rhoSub))
            vAzSub = BzTest*10**-9*np.sqrt((1-alphaRef)/(mu0*rhoSub))
            
            vAxSub = vAxSub/1000
            vAySub = vAySub/1000
            vAzSub = vAzSub/1000
            
            return vAxTest,vAyTest,vAzTest,vAxSub,vAySub,vAzSub
        else:        
            return vAxTest,vAyTest,vAzTest

            
    def WalenPlot(self,fig,ax,mode,RHSx,RHSy,RHSz,LHSx,LHSy,LHSz):
        #Plot Walen Relation for input LHS/RHS for all 3-components at once
        #Run once figure has been constructed
        
        #Define figure properties for light and dark mode
        if mode == 'light':
            ctick = 'k' #Black ticks/labels
        elif mode == 'dark':
            ctick = 'w' #White ticks/labels
        else: #raise exception for incorrect input
            raise Exception('Incorrect Colour Mode')
        
        #Plot x and y axis parameters
        ax.scatter(RHSx,LHSx,marker='.',label='x') #use different markers for each component
        ax.scatter(RHSy,LHSy,marker='+',label='y')
        ax.scatter(RHSz,LHSz,marker='x',label='z')
        ax.legend(borderpad=0.5,labelspacing=0.8,prop={'size': 16,'family':'sans-serif'},edgecolor=ctick)

        #Add general linear gradient to plotted data
        def line_func(x,a,b):
            y = a*x+b
            return y
        
        #Optimise fit between gradient and data 
        fit = optimize.curve_fit(line_func,xdata=np.concatenate((RHSx,RHSy,RHSz)),
                                 ydata=np.concatenate((LHSx,LHSy,LHSz)))
        
        a = fit[0][0] #optimised gradient
        b = fit[0][1] #optimised y-intercept
        pcov = fit[1] #optimised 2x2 covariance matrix
        
        #Calculate correlation coefficient
        corr = np.corrcoef(np.concatenate((RHSx,RHSy,RHSz)),np.concatenate((LHSx,LHSy,LHSz)))[0,1]
        
        #Set axes limit for each plot
        hor_max = np.concatenate((RHSx,RHSy,RHSz)).max()
        hor_min = np.concatenate((RHSx,RHSy,RHSz)).min()
        ver_max = np.concatenate((LHSx,LHSy,LHSz)).max()
        ver_min = np.concatenate((LHSx,LHSy,LHSz)).min()
        lim_min = np.min((hor_min,ver_min))
        lim_max = np.max((hor_max,ver_max))
        
        #Set plot limts
        ax.set_xlim([lim_min,lim_max])
        ax.set_ylim([lim_min,lim_max])
        
        #Create y=x lines across the plot
        line_arr = np.linspace(lim_min,lim_max,10) #set 10 points for plotting
        if corr >= 0:
            ax.plot(line_arr,line_arr,color=ctick,linewidth='1')
        elif corr < 0:
            ax.plot(line_arr,-1*line_arr,color=ctick,linewidth='1')
        
        
        return fig,a,b,corr
    
    def csv_recorder2(self):
        #Method to record plotted data for 2-panels Walen test
        
        start_dtstr = self.test_dti.strftime('%Y-%m-%dT%H-%M-%SZ')
        stop_dtstr = self.test_dtf.strftime('%Y-%m-%dT%H-%M-%SZ')
        csv_filename = '_'.join(['Walen',start_dtstr,stop_dtstr]) + '.csv'
        with open(csv_filename,'w', newline = '') as csvfile:
            csvwriter = csv.writer(csvfile)
            for i in np.arange(len(self.FGMTestdatetime)):
                csvwriter.writerow([self.FGMTestdatetime[i],
                                    self.EHTxTest[i],self.EHTyTest[i],self.EHTzTest[i],
                                    self.ECxTest[i],self.ECyTest[i],self.ECzTest[i],
                                    self.vAxTest[i],self.vAyTest[i],self.vAzTest[i],
                                    self.vdiffxTest[i],self.vdiffyTest[i],self.vdiffzTest[i]])
        csvfile.close()
    
    def csv_recorder3(self):
        #Method to record plotted data for 3-panels Walen test
        
        start_dtstr = self.test_dti.strftime('%Y-%m-%dT%H-%M-%SZ')
        stop_dtstr = self.test_dtf.strftime('%Y-%m-%dT%H-%M-%SZ')
        csv_filename = '_'.join(['Walen',start_dtstr,stop_dtstr]) + '.csv'
        with open(csv_filename,'w', newline = '') as csvfile:
            csvwriter = csv.writer(csvfile)
            for i in np.arange(len(self.FGMTestdatetime)):
                csvwriter.writerow([self.FGMTestdatetime[i],
                                    self.EHTxTest[i],self.EHTyTest[i],self.EHTzTest[i],
                                    self.ECxTest[i],self.ECyTest[i],self.ECzTest[i],
                                    self.vAxTest[i],self.vAyTest[i],self.vAzTest[i],
                                    self.vdiffxTest[i],self.vdiffyTest[i],self.vdiffzTest[i],
                                    self.vAxSub[i],self.vAySub[i],self.vAzSub[i]])
            csvfile.close()
        
          
        
"""
Script for Application
"""

test_dti = dt.datetime(2002,3,18,14,57,25)
test_dtf = dt.datetime(2002,3,18,15,3,6)
ref_dti = dt.datetime(2002,3,18,15,15,0)
ref_dtf = dt.datetime(2002,3,18,15,17,0) #need interval long enough to have 2+ data points. Code only take first data point in this range.
SC = 'C1'

#with density substution
W = Walen(SC,test_dti,test_dtf,ref_dti,ref_dtf,sub_only=True,title=False,mode='dark',saveformat='png',data_record=True) 
#without density substition
#W = Walen(SC,test_dti,test_dtf)

