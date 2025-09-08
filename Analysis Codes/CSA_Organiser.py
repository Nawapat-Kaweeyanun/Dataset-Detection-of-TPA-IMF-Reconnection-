"""
Edition Date: 2024-August-08
@author: Nawapat Kaweeyanun
"""

"""
Objective: Reorganise Cluster Science Archives downloaded folders into year-month-hierarchy
"""

import os
import shutil


"""
Recommend: Move all downloaded CSA folders (NOT CDF files) into a bigger directory
"""

#folder where all CSA folders are gathered
sourcefolder = 'CSA_ToMove'

#folder where all reorganised files will be deposited
destifolder = 'Datafiles_ClusterFGM'


def fileshifter(sourcefolder,destifolder):

    #create list of all CSA folders
    CSAFolderlist = os.listdir(sourcefolder)
    CSAFolderlist.sort()

    for CSAFolder in CSAFolderlist:
        CSAPath = sourcefolder + '/' + CSAFolder
        DataTypeList = os.listdir(CSAPath)
        DataTypeList.sort()

        for DataType in DataTypeList:
            DataPath = CSAPath + '/' + DataType
            Filelist = os.listdir(DataPath)
            Filelist.sort()
        
            if len(Filelist) == 0: #no files
                continue
    
            for File in Filelist:
                #construct source path
                sourcepath = '{}/{}/{}/{}'.format(sourcefolder,CSAFolder,DataType,File)
            
                #determine path of destination DIRECTORY
                #print(File)
                date = File.split('__')[1].split('_')[0] #use start date
                yr = int(date[0:4])
                mo = int(date[4:6])
                destidir = '{}/{}/{}/{}/'.format(destifolder,DataType,yr,mo)
                os.makedirs(destidir,exist_ok=True)
                shutil.move(sourcepath,destidir) #move into directory
                