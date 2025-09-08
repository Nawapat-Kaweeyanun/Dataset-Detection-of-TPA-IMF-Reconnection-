"""
Edition Date: 2024-December-19
@author: Nawapat Kaweeyanun
"""

"""
Objective: Define function for download and unpacking Cluster Science Archive data, to be called for each file

"""

import csv
import os
from requests import get
import tarfile
import numpy as np
import shutil
import CSA_Organiser as Corg

"""
Use request to download file.
url: 'https://csa.esac.esa.int/csa-sl-tap/data'
params = {'RETRIEVAL_TYPE': 'product','DATASET_ID': DataID,
          'START_DATE': StartTime,'END_DATE': EndTime,
          'DELIVERY_FORMAT': 'CDF','DELIVERY_INTERVAL': 'hourly'}
See PanelPlotter for Application
"""
 
def download(url, params, file_name):
    with open(file_name, "wb") as file: #open in binary mode
        response = get(url, params=params) #get request
        file.write(response.content) #write to file

"""
#Unpack tar file, obtain CDF under path: downloadname/DataID/filename
"""
def unpack_tgz(file_name):
    with tarfile.open(file_name) as tar:
        tarname = tar.getnames()
        tar.extractall()
        cdffilepath = tarname[0]
    return cdffilepath    
    

"""
Example script for downloading IMAGE data within specified time period

(Un-comment to use)
"""
"""

#Set download period datetime
StartTimes = ['2002-03-18T14:00:00Z']
EndTimes = ['2002-03-18T15:30:00Z']

#Set target URL, target data ID, and directory to which files are being downloaded
myurl = 'https://csa.esac.esa.int/csa-sl-tap/data'
DataID = 'CC_CP_AUX_ECLAT_IMAGE_WIC' #fi ===for DataSetID category 
MoveTarget = 'Datafiles_ClusterFGM'  

#Loop over all specified periods
for n in np.arange(len(StartTimes)): #len(StartTimes)
    
    #Set download parameter dictionary
    query_specs = {'RETRIEVAL_TYPE': 'product',
               'DATASET_ID': DataID,
               'START_DATE': StartTimes[n],
               'END_DATE': EndTimes[n],
               'DELIVERY_FORMAT': 'CDF',
               'DELIVERY_INTERVAL': 'hourly'}

    #Activate download function
    download(myurl, query_specs, 'tap_download.tgz')

    #Untar & unzip the file, then obtain file name
    cdffilepath = unpack_tgz('tap_download.tgz')
    print('Start time ' + str(n+1) + ' Downloaded') 
    
    #Create intermediate data directory. Remove previous run if exists.
    if 'CSA_ToMove' in os.listdir(os.getcwd()):
        shutil.rmtree('CSA_ToMove')
        os.makedirs('CSA_ToMove', exist_ok=True)
    
    #Move folder to intermediate folder
    filedir = '/'.join(cdffilepath.split('/')[0:-1])
    shutil.move(filedir,'CSA_ToMove' + '/' + filedir)
    
#Remove CSA Download Folders to reduce clutter
DownloadFolderList = [d for d in os.listdir(os.getcwd()) if 'CSA_Download' in d]
for DLF in DownloadFolderList:
    try:
        shutil.rmtree(DLF) #remove dir and files inside
    except: #skip a few parsepoint error files.
        continue
        
#Shift files to final directory destination
Corg.fileshifter('CSA_ToMove',MoveTarget)
"""