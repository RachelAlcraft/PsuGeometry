

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/1000Structures/'


# hand crafted report for 2bw4 non planar omega
from PsuGeometry import GeoReport as geor
import pandas as pd


##############################################################################################
### USER CHOICE loading or saing data and whether to override #################################
###############################################################################################
mode = 'ALL' # LOAD, SAVE or ALL
override = False
ghost = False
loadDssp = False
loadEd = True
###############################################################################################


pdbList = []
pdbListPath = printPath + '1000 Structures.csv'
pdbdata = pd.read_csv(pdbListPath)
pdbList = pdbdata['pdb_code']


if mode=='SAVE': #My computer cannot cope with the memory requirements!!

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
            #specify the paramaters
            geoList = ['N:CA:C','PSI']
            hueList = ['dssp','aa','bfactor','2FoFc','rid','resolution'] # note the hues are the sum of the atoms
            #Create the dataframe
            data = georep.getGeoemtryCsv(geoList, hueList)
            #And save this to file
            data.to_csv(fullFileName, index=False)


elif mode=='LOAD': # then load the csv files up to produce the reports
    isEmpty = True
    allData = None
    # Create an empty report object
    georep = geor.GeoReport([],pdbDataPath,edDataPath,printPath,ed=loadEd,dssp=loadDssp)
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
    georep.addScatter(data=allData, geoX='aa', geoY='bfactor', title='Info', hue='rid',palette='rainbow')
    georep.addScatter(data=allData, geoX='resolution', geoY='aa', title='Info', hue='resolution')

    georep.addScatter(data=allData, geoX='N:CA:C', geoY='N:CA:C:N+1',title='TAU-PSI',  hue='resolution')
    georep.addProbability(data=allData, geoX='N:CA:C', geoY='N:CA:C:N+1', title='TAU-PSI',palette='cubehelix_r')
    georep.addHistogram(data=allData, geoX='N:CA:C', title='TAU')

    # Print the report
    georep.printToHtml('1000 High Res Structures', 3, '1000')

else:
    georep = geor.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=True)

    geoList = ['TAU', 'PSI','N:O','CB:O','PHI']
    hueList = ['aa', 'rid', 'resolution','dssp']  # note the hues are the sum of the atoms
    # Create the dataframe
    data = georep.getGeoemtryCsv(geoList, hueList)

    print('Creating report')
    georep.addScatter(data=data, geoX='PSI', geoY='N:O', title='PSI-N:O', hue='resolution')
    georep.addScatter(data=data, geoX='PSI', geoY='CB:O', title='PSI-CB:O', hue='resolution')
    georep.addScatter(data=data, geoX='N:O', geoY='CB:O', title='N:O-CB:O', hue='resolution')

    georep.addScatter(data=data, geoX='PSI', geoY='N:O', title='PSI-N:O', hue='dssp',palette='tab20')
    georep.addScatter(data=data, geoX='PSI', geoY='CB:O', title='PSI-CB:O', hue='dssp',palette='tab20')
    georep.addScatter(data=data, geoX='N:O', geoY='CB:O', title='N:O-CB:O', hue='dssp',palette='tab20')

    georep.addScatter(data=data, geoX='PHI', geoY='PSI', title='Ramachandran Plot', hue='dssp',palette='tab20')
    georep.addScatter(data=data, geoX='resolution', geoY='pdbCode', title='Info', hue='resolution')
    georep.addHistogram(data=data, geoX='N:CA:C', title='TAU')

    georep.addScatter(data=data,geoX='N:CA:C', geoY='N:CA:C:N+1', title='TAU-PSI', hue='resolution')
    georep.addProbability(data=data,geoX='N:CA:C', geoY='N:CA:C:N+1', title='TAU-PSI', palette='cubehelix_r')
    georep.addScatter(data=data, geoX='N:CA:C', geoY='N:CA:C:N+1', title='TAU-PSI', hue='dssp',palette='tab20')

    georep.addDifference(dataA=data, dataB=data, geoX='N:CA:C', geoY='N:CA:C:N+1', restrictionsA={'aa': 'VAL'},restrictionsB={'aa': 'THR'})






    # Print the report
    georep.printToHtml('1000 High Res Structures', 3, '1000')