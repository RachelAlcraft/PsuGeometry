

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/bimodal/'


# hand crafted report for 2bw4 non planar omega
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop
import pandas as pd
from PsuGeometry import GeoPlot as geopl

##############################################################################################
### USER CHOICE loading or saing data and whethr to override #################################
###############################################################################################
mode = 'LOAD' # LOAD, SAVE or ALL
override = False
ghost = False
loadDssp = False
loadEd = False
###############################################################################################


pdbList = []
pdbListPath = printPath + 'pdbs.csv'
pdbdata = pd.read_csv(pdbListPath)
pdbList = pdbdata['pdb_code']


if mode=='SAVE': #My computer cannot cope with the memory requirements!!

    count = 0
    num = len(pdbList)
    for pdb in pdbList:
        count += 1
        pdb = pdb.lower()
        fileName = pdb + '_tau.csv'
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
        fileName = pdb + '_tau.csv'
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

    for aa in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']:
        print('Creating',aa)
        georep.addScatter(data=allData, geoX='N:CA:C', geoY='N:CA:C:N+1', title='', hue='resolution',restrictions={'aa':aa},ghost=ghost)
        georep.addProbability(data=allData, geoX='N:CA:C', geoY='N:CA:C:N+1', title='', palette='cubehelix_r',restrictions={'aa':aa},ghost=ghost)
        georep.addHistogram(data=allData, geoX='N:CA:C', title='',restrictions={'aa':aa},ghost=ghost)
        georep.addDifference(dataA=allData, dataB=allData, geoX='N:CA:C', geoY='N:CA:C:N+1', restrictionsA={'aa': aa}, exclusionsB={'aa': aa})

    # Print the report
    georep.printToHtml('Bimodal Tau', 3, 'tau')

else:
    georep = geor.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=loadEd, dssp=loadDssp)

    geoList = ['TAU', 'PSI']
    hueList = ['aa', 'rid', 'resolution']  # note the hues are the sum of the atoms
    # Create the dataframe
    data = georep.getGeoemtryCsv(geoList, hueList)

    print('Creating report')
    georep.addScatter(data=data,geoX='resolution', geoY='pdbCode', title='Info', hue='resolution')
    georep.addScatter(data=data,geoX='rid', geoY='pdbCode', title='Info', hue='resolution')
    georep.addScatter(data=data,geoX='resolution', geoY='aa', title='Info', hue='resolution')

    georep.addScatter(data=data,geoX='N:CA:C', geoY='N:CA:C:N+1', title='TAU-PSI', hue='resolution')
    georep.addProbability(data=data,geoX='N:CA:C', geoY='N:CA:C:N+1', title='TAU-PSI', palette='cubehelix_r')
    georep.addHistogram(data=data,geoX='N:CA:C', title='TAU')

    for aa in ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']:
        georep.addScatter(data=data,geoX='N:CA:C', geoY='N:CA:C:N+1', title='', hue='resolution',restrictions={'aa': aa}, ghost=ghost)
        georep.addProbability(data=data,geoX='N:CA:C', geoY='N:CA:C:N+1', title='', palette='cubehelix_r',restrictions={'aa': aa}, ghost=ghost)
        georep.addHistogram(data=data,geoX='N:CA:C', title='', restrictions={'aa': aa}, ghost=ghost)
        georep.addDifference(data=data,geoX='N:CA:C', geoY='N:CA:C:N+1', restrictionsA={'aa': aa},exclusionsB={'aa': aa})

    # Print the report
    georep.printToHtml('Bimodal Tau', 3, 'tau')