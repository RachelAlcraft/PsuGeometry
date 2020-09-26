
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/jask/'


# hand crafted report for 2bw4 non planar omega
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop
import pandas as pd
from PsuGeometry import GeoPlot as geopl

pdbList = ['1ejg', '1ucs', '1us0', '1yk4', '1r6j', '1hje', '3al1', '2b97', '1gci', '1x6z']
saveToCsv = False
override = False

if saveToCsv: #Jaskolski test My computer cannot cope with the memory requirements!!

    for pdb in pdbList:
        pdb = pdb.lower()
        fileName = pdb + '_jask.csv'
        fullFileName = printPath + fileName
        needToCreate = override
        data = None
        try:
            data = pd.read_csv(fullFileName)
        except:
            print('File does not exist, creating...', fullFileName)
            needToCreate = True

        if needToCreate:
            #Init the report object
            georep = geor.GeoReport([pdb],pdbDataPath,edDataPath,printPath)
            #specify the paramaters
            geoList = ['C-1:N','N:CA','CA:C','C:N+1','C:O','N:CA:C','N:CA:C:N+1']
            hueList = ['dssp','aa','bfactor','2FoFc','rid','resolution'] # note the hues are the sum od the atoms
            #Create the dataframe
            data = georep.getGeoemtryCsv(geoList, hueList)
            #And save this to file
            data.to_csv(fullFileName, index=False)


else: # then load the csv files up to produce the reports
    isEmpty = True
    allData = None
    # Create an empty report object
    georep = geor.GeoReport([],pdbDataPath,edDataPath,printPath)
    for pdb in pdbList:
        pdb = pdb.lower()
        fileName = pdb + '_jask.csv'
        fullFileName = printPath + fileName
        needToCreate = override
        data = None
        try:
            data = pd.read_csv(fullFileName)
            if isEmpty:
                allData = data
                isEmpty = False
            else:
                allData = pd.concat([allData,data])
        except:
            print('File does not exist',fullFileName)

    allData = allData.query('aa !="ACE"')
    allDataexcGlyPro = allData.query('aa !="PRO" and aa!="GLY"') #N:CA
    allDataexcGly = allData.query('aa!="GLY"') #CA:C
    allDataexcPro = allData.query('aa!="PRO"')  # C-1:N
    allDataGly = allData.query('aa=="GLY"') #TAU
    allDataPro = allData.query('aa =="PRO"') # TAU
    printList = []
    georep.addScatter(data=allData, geoX='pdbCode', geoY='resolution',title='Resolutions',  hue='pdbCode', palette='Accent')
    #georep.addHistogram(data=allData, geoX='pdbCode',title='PDBs')
    #georep.addScatter(data=allData, geoX='pdbCode', geoY='aa', title='Amino Acids', hue='aa', palette='Accent')
    georep.addProbability(data=allData, geoX='N:CA:C', geoY='N:CA:C:N+1', title='TAU-PSI',palette='cubehelix_r')
    georep.addHistogram(data=allDataexcPro, geoX='C-1:N',title='C-1:N exc PRO')
    georep.addHistogram(data=allDataexcGlyPro, geoX='N:CA',title='N:CA exc GLY PRO')
    georep.addHistogram(data=allDataexcGly, geoX='CA:C',title='CA:C exc GLY')
    georep.addHistogram(data=allData, geoX='C:O',title='C:O')
    georep.addHistogram(data=allDataexcGlyPro, geoX='N:CA:C',title='N:CA:C exc GLY PRO')
    georep.addHistogram(data=allDataGly, geoX='N:CA:C',title='N:CA:C GLY')
    georep.addHistogram(data=allDataPro, geoX='N:CA:C',title='N:CA:C PRO')
    #georep.addHistogram(data=allData,geoX='N:CA:C',splitKey='aa',title='TAU')
    #georep.addHistogram(data=allDataexcPro, geoX='C-1:N',title='C-1:N exc PRO', splitKey='pdbCode')
    #georep.addHistogram(data=allDataexcGlyPro, geoX='N:CA',title='N:CA exc GLY PRO', splitKey='pdbCode')
    #georep.addHistogram(data=allDataexcGly, geoX='CA:C',title='CA:C exc GLY', splitKey='pdbCode')
    #georep.addHistogram(data=allData, geoX='C:O',title='C:O', splitKey='pdbCode')
    # Print the report
    georep.printToHtml('Jaskolski Pdbs', 3, 'jaskolski')
