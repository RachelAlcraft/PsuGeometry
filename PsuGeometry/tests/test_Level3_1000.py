# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdbLists as geol
import pandas as pd
'''
This script produces correlation reports for a high resolution dataset of 1000 structures.
Due to memory constraints in loading this large amount of data in python (well for me on my cheap laptop)
the script is set up to output csv files on a per structure basis as csv files in the first pass.
This means it can keep restarting every time it runs out of memory :-)
This is not the way we write code :-) 

The potential reason for memory use is the loading of the electron density, and the use of electron density as a hue.
Without that, it would not use all the memory. 
In another language, we could delete the electron density data as we go! Somebody help me do this with python...

This script dies not use the electron density because it would take for ever, see Level3_test_Electron|Density
But the pattern is set up for slow or large amounts of data.
Once seriablied to CSV it is quite quick to run reports.

First, keep running it in mode SAVE until it is finished.
Then run it in mode LOAD, loading simple csv files means the data can easily be processed, even on that large scale, 
producing an elegant, data-rich report.
'''

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/1000Structures/'

##############################################################################################
### USER CHOICE loading or saing data and whether to override #################################
###############################################################################################
mode = 'LOAD' # LOAD or SAVE
override = True # choose if you want to write over the already saved files, if not it will start from where it left off each time
ghost = False # choose if you want to have the silver ghost structure in the report
loadDssp = True # whether or not to load the secondary structure library - some memory implications
loadEd = False # whether or not to load the electron density, it will be much more efficient of you do not (but you won't have the electron density!). Will use bfactor instead of ed if you turn this off.
geoList = ['TAU','PSI','PHI','OMEGA','N:O','CB:O','N:CA','CA:C','C:O','C:N+1','CA:C:N+1','N+1:C:O','CA:C:O']
hueList = ['dssp','aa','bfactor','2FoFc','rid','resolution']
###############################################################################################

pdbList = []
pdbList = geol.GeoPdbLists().getList1000()

if mode=='SAVE': # Then we are going to write the csv files out, and keep re-starting until we are done
    count = 0
    num = len(pdbList)
    for pdb in pdbList:
        count += 1
        pdb = pdb.lower()
        fileName = pdb + '_1000.csv'
        fullFileName = printPath + fileName
        needToCreate = override
        data = None
        try:
            data = pd.read_csv(fullFileName)
            print('Exists', fullFileName)
        except:
            print('File does not exist, creating...', fullFileName)
            needToCreate = True

        if needToCreate:
            #Init the report object
            print('Loading',count,'/',num)
            georep = geor.GeoReport([pdb],pdbDataPath,edDataPath,printPath,ed=loadEd,dssp=loadDssp)
            #Create the dataframe
            data = georep.getGeoemtryCsv(geoList, hueList)
            #And save this to file
            data.to_csv(fullFileName, index=False)


elif mode=='LOAD': # then load the csv files up to produce the reports
    isEmpty = True
    allData = None
    # Create an empty report object
    georep = geor.GeoReport([],pdbDataPath,edDataPath,printPath,includePdbs=False,ed=loadEd,dssp=loadDssp)
    count = 0
    dfs = []
    for pdb in pdbList:
        pdb = pdb.lower()
        count = count + 1
        fileName = pdb + '_1000.csv'
        fullFileName = printPath + fileName
        needToCreate = override

        data = None
        try:
            print('Loading', count, '/', len(pdbList))
            data = pd.read_csv(fullFileName)
            dfs.append(data)
        except:
            print('File does not exist',fullFileName)

    print('Concatenating')
    allData = pd.concat(dfs, ignore_index=True)
    print('Creating report')
    georep.addScatter(data=allData, geoX='resolution', geoY='pdbCode', title='Info', hue='resolution')
    georep.addScatter(data=allData, geoX='PHI', geoY='PSI', title='Ramachandran', hue='bfactor',palette='ocean_r',sort='ASC')
    georep.addScatter(data=allData, geoX='OMEGA', geoY='TAU', title='Omega-Tau', hue='dssp',palette='nipy_spectral',sort='NON')

    georep.addScatter(data=allData, geoX='PHI', geoY='PSI', title='Ramachandran', hue='dssp',palette='nipy_spectral',sort='NON')
    georep.addScatter(data=allData, geoX='TAU', geoY='PSI', title='Tau-Psi', hue='dssp',palette='nipy_spectral',sort='NON')
    georep.addScatter(data=allData, geoX='PHI', geoY='TAU', title='Phi-Tau', hue='dssp', palette='nipy_spectral',sort='NON')

    georep.addScatter(data=allData, geoX='PSI', geoY='N:O', title='Psi-N:O', hue='dssp', palette='nipy_spectral',sort='NON')
    georep.addScatter(data=allData, geoX='PSI', geoY='CB:O', title='Psi-CB:O', hue='dssp', palette='nipy_spectral',sort='NON')
    georep.addScatter(data=allData, geoX='N:O', geoY='CB:O', title='N:O-CB:O', hue='dssp', palette='nipy_spectral',sort='NON')

    georep.addScatter(data=allData, geoX='N:CA:C', geoY='N:CA:C:N+1',title='TAU-PSI',  hue='resolution')
    georep.addProbability(data=allData, geoX='N:CA:C', geoY='N:CA:C:N+1', title='TAU-PSI',palette='cubehelix_r')
    georep.addHistogram(data=allData, geoX='N:CA:C', title='TAU')

    georep.addHistogram(data=allData, geoX='CA:C', title='CA:C',exclusions={'pdbCode':'4y9w'})
    georep.addHistogram(data=allData, geoX='C:O', title='C:O',exclusions={'pdbCode':'1d5t'})
    georep.addHistogram(data=allData, geoX='C:N+1', title='C:N+1',hue='rid',exclusions={'pdbCode':'3aj4','pdbCode':'2cnu'})

    georep.addHistogram(data=allData, geoX='CA:C:N+1', title='CA:C:N+1')
    georep.addHistogram(data=allData, geoX='N+1:C:O', title='N:C:O')
    georep.addHistogram(data=allData, geoX='CA:C:O', title='CA:C:O')

    # Print the report
    georep.printToHtml('1000 High Res Structures', 3, '1000')

