"""
Edition Date: 2025-September-01
@author: Nawapat Kaweeyanun
"""

"""
Objective: Define functions to search for correct cdf IMAGE file for a given datetime and camera

Prerequisite: 
    IMAGE files must be saved in FolderName/DataID/year/month structure.
    SSUSI files must be saved in FolderName/dataN/Spacecraft/apl/edr-aur/year/day-of-year structure.
    GUVI files must be saved in FolderName/edr-aur/year/day-of-year structure.
    The SSUSI_GUVI_Data extractor must be in the same directory as this file.
"""

import os
import datetime as dt
import SSUSI_GUVI_Data as SGD



def IMAGE_filesearch(IMFolder,IMtype,dtime): #search function for IMAGE
    
    #Determine data ID for WIC and S12 (SI-12)
    if IMtype == 'WIC':
        DataID = 'CC_CP_AUX_ECLAT_IMAGE_WIC'
    elif IMtype == 'S12':
        DataID = 'CC_CP_AUX_ECLAT_IMAGE_S12'
        
    #Obtain directory path
    destifolder = '{}/{}/{}/{}'.format(IMFolder,DataID,dtime.year,dtime.month)
    
    try:
        #Sort files in destination folder alphabetically
        filelist = os.listdir(destifolder)
        filelist.sort()
        
        #Construct date and time strings from given datetime
        date_str = str(dtime).split(' ')[0].replace('-','')
        time_str = str(dtime).split(' ')[1].replace(':','')
        
        #Pick file that match start date with correct time range
        for filename in filelist:
            startdate_str = filename.split('__')[1].split('_')[0]
            start_hour = int(filename.split('__')[1].split('_')[1][0:2])
            end_hour = int(filename.split('__')[1].split('_')[3][0:2])
            if date_str == startdate_str: #right date
                if dtime.hour < 23:
                    if start_hour <= dtime.hour and end_hour > dtime.hour:
                        target_file = filename
                    else:
                        continue
                elif dtime.hour == 23:
                    if start_hour == 23:
                        target_file = filename
                    else:
                        continue
            
            #check if correct file is found. If so, construct path to that file
            if 'target_file' not in locals():
                target_filepath = None
            else:
                target_filepath = destifolder + '/' + target_file
    except:
        print('File list not generated for folder {} - check if data exist.'.format(destifolder))
        target_filepath = None
                
    return target_filepath

def SSUSI_filesearch(SSFolder,SStype,SSinstr,dtime):
    #Function to search for a SSUSI scan from orbit closest to input datetime.
    
    #Obtain day-of-year from datetime and format it as a string.
    yday = dtime.timetuple().tm_yday
    yday_str = '{:03d}'.format(yday)
    
    #Construct the data file directory path.
    destifolder = '{}/dataN/{}/apl/{}/{}/{}'.format(SSFolder,SSinstr.lower(),
                                                SStype.lower(),dtime.year,
                                                yday_str)
    
    #Use try-except to skip folders that do not contain data.
    try:
        filelist = os.listdir(destifolder)
        filelist.sort()
        
        #Construct date string from datetime.
        date_str = str(dtime).split(' ')[0].replace('-','')

        #Search for file whose time range contains the input datetime.
        for filename in filelist:
            startdate_str = filename.split('_')[6].split('.')[1]
        
            #if dates match, extract file to determine time match
            if date_str == startdate_str:
                filepath = destifolder + '/' + filename
                DataSet = SGD.SGExtract(filepath)
            
                if DataSet.start_dt < dtime and DataSet.stop_dt > dtime:
                    target_file = filename
                
            #if time match, construct path to that file. Otherwise continue and target filepath remains None.
            if 'target_file' not in locals():
                target_filepath = None
            else:
                #construct path to that file
                target_filepath = destifolder + '/' + target_file
                    
    #If there is no match at all, print exception and check if data exist.
    except Exception:
        print('File list not generated for folder {} - check if data exist.'.format(destifolder))  
        target_filepath = None
   
    return target_filepath


def GUVI_filesearch(GUFolder,GUtype,dtime):
    #Function to search for a GUVI scan from orbit closest to input datetime.
    
    #Obtain day-of-year
    yday = dtime.timetuple().tm_yday
    yday_str = '{:03d}'.format(yday)
        
    #Construct the data file directory path.
    destifolder = '{}/{}/{}/{}'.format(GUFolder,GUtype.lower(),
                                          dtime.year,yday_str)
    
    #Use try-except to skip folders that do not contain data.
    try:
        filelist = os.listdir(destifolder)
        filelist.sort()
        
        for filename in filelist:
            #Convert start/end datetime strings of each file to datetime objects.
            #Because datetime is already in filename, no need to do extraction like SSUSI.
            start_dtstr = filename.split('_')[3].split('-')[0]
            start_dt = dt.datetime.strptime(start_dtstr,'%Y%j%H%M%S')
            stop_dtstr = filename.split('_')[3].split('-')[1]
            stop_dt = dt.datetime.strptime(stop_dtstr,'%Y%j%H%M%S')
            
            #If a file's datetime range matches input datetime, construct a filepath.
            if dtime > start_dt and dtime < stop_dt:
                target_filepath = destifolder + '/' + filename
                break
            else:
                target_filepath = None#
                continue
    #If there is no match at all, print exception and check if data exist.
    except:
        print('File list not generated for folder {} - check if data exist.'.format(destifolder))  
        target_filepath = None
   
    return target_filepath
    