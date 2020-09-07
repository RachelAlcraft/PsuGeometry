
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/Jask/'


# hand crafted report for 2bw4 non planar omega
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop
import pandas as pd

pdbList = ['1ejg', '1ucs', '1us0', '1yk4', '1r6j', '1hje', '3al1', '2b97', '1gci', '1x6z']
saveToCsv = False
override = True

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
            geoPdb = geop.GeoPdb(pdb, pdbDataPath,edDataPath)
            #Init the report object
            georep = geor.GeoReport([geoPdb])
            #specify the paramaters
            geoList = ['C-1:N','N:CA','CA:C','C:N+1','C:O','N:CA:C']
            hueList = ['dssp','aa','bfactor','2FoFc','rid','resolution'] # note the hues are the sum od the atoms
            #Create the dataframe
            data = georep.getGeoemtryCsv(geoList, hueList)
            #And save this to file
            data.to_csv(fullFileName, index=False)


else: # then load the csv files up to produce the reports
    isEmpty = True
    allData = None
    # Create an empty report object
    georep = geor.GeoReport([])
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
    printList.append(['Resolutions', allData, 'pdbCode', 'resolution', '', 'pdbCode', 'Accent', False, 0, 0])
    printList.append(['PDBs', allData, 'pdbCode', ''])
    printList.append(['Amino Acids', allData, 'pdbCode', 'aa', '', 'aa', 'Accent', False, 0, 0])
    #printList.append(['Amino Acids', allData, 'aa', ''])
    printList.append(['C-1:N exc PRO', allDataexcPro, 'C-1:N', ''])
    printList.append(['N:CA exc GLY PRO', allDataexcGlyPro, 'N:CA', ''])
    printList.append(['CA:C exc GLY', allDataexcGly, 'CA:C', ''])
    printList.append(['C:O', allData, 'C:O', ''])
    printList.append(['N:CA:C exc GLY PRO', allDataexcGlyPro, 'N:CA:C', ''])
    printList.append(['N:CA:C GLY', allDataGly, 'N:CA:C', ''])
    printList.append(['N:CA:C PRO', allDataPro, 'N:CA:C', ''])
    printList.append(['TAU', allData, 'N:CA:C', '','aa'])
    printList.append(['C-1:N exc PRO', allDataexcPro, 'C-1:N', '','pdbCode'])
    printList.append(['N:CA exc GLY PRO', allDataexcGlyPro, 'N:CA', '','pdbCode'])
    printList.append(['CA:C exc GLY', allDataexcGly, 'CA:C', '','pdbCode'])
    printList.append(['C:O', allData, 'C:O', '','pdbCode'])
    # Print the report
    georep.printCsvToHtml(printList, [], 'Jaskolski Pdbs', 3, printPath, 'jaskolski')