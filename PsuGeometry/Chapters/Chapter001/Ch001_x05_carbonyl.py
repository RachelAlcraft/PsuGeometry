'''
This script creates a file for hand chosen geos, so it is slower as they are no serialised.
This particular report is designed to look at hydrogen bonding on the carbonyl oxygen
'''

import pandas as pd
from PsuGeometry import GeoReport as psu
import Ch000_Functions as help
import matplotlib
print(matplotlib.__version__)

filesPDBRoot ='C:/Dev/Github/ProteinDataFiles/pdb_data/'
filesADJRoot ='C:/Dev/Github/ProteinDataFiles/pdb_out/Fov2_ADJ/' #adjusted on Fo at 3 degrees thevenaz
loadPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/'
printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'

geos = ['TAU','TAU+1','TAU-1','CA:C:O','O:C:N+1','CA:C:N+1','CA-1:CA:CA+1',
        'O-1:C-1','C-1:N','N:CA','CA:C','C:O','C:N+1','N+1:CA+1','CA+1:C+1','C+1:O+1',
        'PHI','PSI','OMEGA','CA-1:C-1:N:CA',
        'CA-1:CA','CA:CA+1','C-1:C','C:C+1','N-1:N','N:N+1',
        'CA-1:N','CA-1:O-1','O-1:N','C-1:CA','N:C','CA:O','CA:N+1','O:N+1','C:CA+1','N+1:C+1',
        'O-1:CA','N:O','O:CA+1','N+1:O+1','N-1:O-1']


title='Backbone Report'
fileName = 'backbone'

createOrLoad = "CREATE" # CREATE or LOAD
if createOrLoad == "LOAD":
    print('### CREATING csv files ###')
    pdbdata = pd.read_csv('../../PdbLists/Pdbs_70.csv')
    pdbListA = pdbdata['PDB'].tolist()[0:]
    pdbListIn = []
    for pdb in pdbListA:
        import os.path
        if os.path.isfile((filesADJRoot + 'pdb' + pdb + '.ent').lower()):
            pdbListIn.append(pdb.lower())
        else:
            print('No file:',(filesADJRoot + 'pdb' + pdb + '.ent').lower())
    print(pdbListIn)
    print("---- Getting bad atom list--------")
    badAtoms = help.getBadAtomsListFromFile(loadPath, "badatoms.csv")  # Get the bad atoms list we will use to reduce the list further
    print("---- Making unrestricted--------")
    dataPdbUn = help.makeCsv('PDB', pdbListIn, geos, [],True)
    print("---- Making unrestricted--------")
    dataPdbRes = help.makeCsv('PDB', pdbListIn, geos, [],False)
    dataPdbRes = help.applyRestrictions(dataPdbRes)
    print("---- Making reduced--------")
    dataPdbCut = help.makeCsv('PDB', pdbListIn, geos, badAtoms,False)
    dataPdbCut = help.applyRestrictions(dataPdbCut)
    print("---- Making adjusted--------")
    dataPdbAdj = help.makeCsv('ADJUSTED', pdbListIn, geos, badAtoms,False)
    dataPdbAdj = help.applyRestrictions(dataPdbAdj)
    # embellish with dssp - the dssp file was created ages ago from the linux laptop
    pdbdssp = pd.read_csv('C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/CsvGeos_BEST_Set0DSSPALL.csv')
    try:
        pdbdssp['rid'] = pdbdssp['rid'].astype(str)
        pdbdssp['DSSPID'] = pdbdssp['pdbCode'] + pdbdssp['chain'] + pdbdssp['rid']
        pdbdssp = pdbdssp[['DSSPID', 'dssp']]

    except:
        print('empty csv')

    allList = []
    allList.append([dataPdbUn, loadPath + "bb_unrestricted.csv"])
    allList.append([dataPdbRes, loadPath + "bb_restricted.csv"])
    allList.append([dataPdbCut, loadPath + "bb_reduced.csv"])
    allList.append([dataPdbAdj, loadPath + "bb_adjusted.csv"])

    for dataBlob in allList:
        dataCsv = dataBlob[0]
        dataPath = dataBlob[1]
        try:
            dataCsv['rid'] = dataCsv['rid'].astype(str)
            dataCsv['DSSPID'] = dataCsv['pdbCode'] + dataCsv['chain'] + dataCsv['rid']
        except:
            print('empty csv')

        print('### Writing to file',dataPath,'###')
        dataCsv = dataCsv.set_index('DSSPID').join(pdbdssp.set_index('DSSPID'))
        dataCsv.to_csv(dataPath, index=False)

print('### LOADING csv files ###') # bit rubbish but we didn;t change the object references with dssp
dataPdbUn = pd.read_csv(loadPath + "bb_unrestricted.csv")
dataPdbRes = pd.read_csv(loadPath + "bb_restricted.csv")
dataPdbCut = pd.read_csv(loadPath + "bb_reduced.csv")
dataPdbAdj = pd.read_csv(loadPath + "bb_adjusted.csv")

# Prepare a seperate report for the pdb
onePdbUn = dataPdbUn.query("pdbCode == '4r2x'")
onePdbRes = dataPdbRes.query("pdbCode == '4r2x'")
onePdbCut = dataPdbCut.query("pdbCode == '4r2x'")
onePdbAdj = dataPdbAdj.query("pdbCode == '4r2x'")

#Choose a bfactor to cut at
dataPdbUn = dataPdbUn.query("`C:O_avbfactor` <= 15")
dataPdbRes = dataPdbRes.query("`C:O_avbfactor` <= 15")
dataPdbCut = dataPdbCut.query("`C:O_avbfactor` <= 15")
dataPdbAdj = dataPdbAdj.query("`C:O_avbfactor` <= 15")

dataPdbUn = dataPdbUn.query("pdbCode != '4r2x'")
dataPdbRes = dataPdbRes.query("pdbCode != '4r2x'")
dataPdbCut = dataPdbCut.query("pdbCode != '4r2x'")
dataPdbAdj = dataPdbAdj.query("pdbCode != '4r2x'")


print('### Creating scatter files ###')

geoTriosA = [
            ['PHI', 'PSI', 'dssp',True],
            ['PHI', 'PSI', 'C:O',False],
            ['PHI', 'PSI', 'TAU',False],
            ['PHI', 'PSI', 'TAU-1',False],
            ['PHI', 'PSI', 'TAU+1',False],
            ['PSI', 'N:N+1', 'dssp',True],
            ['PSI', 'N:N+1', 'C:O',False],
            ['PSI', 'N:N+1', 'TAU',False],
            ['PSI', 'N:N+1', 'TAU+1',False],
            ['PSI', 'N:N+1', 'TAU-1',False],
            ['TAU+1', 'TAU', 'C:O',False],
            ['TAU','TAU-1', 'C:O',False],
            ['TAU-1','TAU+1', 'C:O',False],
            ['C:O','TAU+1', 'bfactor',False],
           ]

geoTriosB = [['C:O'],
            ['PHI', 'PSI', 'dssp', True],
            ['PHI', 'PSI', 'bfactor', False],
            ['C:O', 'C:N+1', 'TAU+1', False],
            ['C:O', 'C:N+1', 'O:N+1', False],
            ['CA:C', 'C:O', 'CA:O', False],
            ['CA:C', 'C:N+1', 'CA:N+1', False],
            ['C:O', 'C:N+1', 'TAU+1', False],
            ['C:O', 'C:CA+1', 'TAU+1', False],
            ['C:O', 'C:C+1', 'TAU+1', False],
            ['C:O', 'CA-1:CA:CA+1', 'TAU+1', False],
            ['C:O', 'O:N+1', 'TAU+1', False],
            ['C:O', 'O:CA+1', 'TAU+1', False],
            ['TAU', 'TAU+1', 'C:O', False],
            ['C:O', 'TAU+1', 'C:O_avbfactor', False],
            ]


geoTriosCO = [
['C:O', 'TAU+1', 'dssp', True],
['C:O', 'TAU', 'dssp', True],
['C:O', 'TAU-1', 'dssp', True],
['C:O', 'PSI', 'dssp', True],
['C:O', 'PHI', 'dssp', True],
['C:O', 'OMEGA', 'dssp', True],
['C:O', 'CA-1:C-1:N:CA', 'dssp', True],
['C:O', 'CA-1:CA:CA+1', 'dssp', True],
['C:O','O-1:C-1','dssp', True],
['C:O','C-1:N','dssp', True],
['C:O','N:CA','dssp', True],
['C:O','CA:C','dssp', True],
['C:O','C:O','dssp', True],
['C:O','C:N+1','dssp', True],
['C:O','N+1:CA+1','dssp', True],
['C:O','CA+1:C+1','dssp', True],
['C:O','C+1:O+1','dssp', True],
['C:O','CA-1:CA','dssp', True],
['C:O','CA:CA+1','dssp', True],
['C:O','C-1:C','dssp', True],
['C:O','C:C+1','dssp', True],
['C:O','N-1:N','dssp', True],
['C:O','N:N+1','dssp', True],
['C:O','CA-1:N','dssp', True],
['C:O','CA-1:O-1','dssp', True],
['C:O','O-1:N','dssp', True],
['C:O','C-1:CA','dssp', True],
['C:O','N:C','dssp', True],
['C:O','CA:O','dssp', True],
['C:O','CA:N+1','dssp', True],
['C:O','O:N+1','dssp', True],
['C:O','C:CA+1','dssp', True],
['C:O','N+1:C+1','dssp', True],
['C:O','O-1:CA','dssp', True],
['C:O','N:O','dssp', True],
['C:O','O:CA+1','dssp', True],
['C:O','N+1:O+1','dssp', True],
['C:O','N-1:O-1','dssp', True],
]

geoTriosTAU1 = [
['C:O', 'TAU+1', 'TAU+1', False],
['C:O', 'TAU', 'TAU+1', False],
['C:O', 'TAU-1', 'TAU+1', False],
['C:O', 'PSI', 'TAU+1', False],
['C:O', 'PHI', 'TAU+1', False],
['C:O', 'OMEGA', 'TAU+1', False],
['C:O', 'CA-1:C-1:N:CA', 'TAU+1', False],
['C:O', 'CA-1:CA:CA+1', 'TAU+1', False],
['C:O','O-1:C-1', 'TAU+1', False],
['C:O','C-1:N', 'TAU+1', False],
['C:O','N:CA', 'TAU+1', False],
['C:O','CA:C', 'TAU+1', False],
['C:O','C:O', 'TAU+1', False],
['C:O','C:N+1', 'TAU+1', False],
['C:O','N+1:CA+1', 'TAU+1', False],
['C:O','CA+1:C+1', 'TAU+1', False],
['C:O','C+1:O+1', 'TAU+1', False],
['C:O','CA-1:CA', 'TAU+1', False],
['C:O','CA:CA+1', 'TAU+1', False],
['C:O','C-1:C', 'TAU+1', False],
['C:O','C:C+1', 'TAU+1', False],
['C:O','N-1:N', 'TAU+1', False],
['C:O','N:N+1', 'TAU+1', False],
['C:O','CA-1:N', 'TAU+1', False],
['C:O','CA-1:O-1', 'TAU+1', False],
['C:O','O-1:N', 'TAU+1', False],
['C:O','C-1:CA', 'TAU+1', False],
['C:O','N:C', 'TAU+1', False],
['C:O','CA:O', 'TAU+1', False],
['C:O','CA:N+1', 'TAU+1', False],
['C:O','O:N+1', 'TAU+1', False],
['C:O','C:CA+1','TAU+1', False],
['C:O','N+1:C+1', 'TAU+1', False],
['C:O','O-1:CA', 'TAU+1', False],
['C:O','N:O', 'TAU+1', False],
['C:O','O:CA+1', 'TAU+1', False],
['C:O','N+1:O+1', 'TAU+1', False],
['C:O','N-1:O-1', 'TAU+1', False],
]

'''
help.trioReports(["Unrestricted",dataPdbUn],
                 ["Restricted",dataPdbRes],
                 ["Reduced",dataPdbCut],
                 ["Adjusted",dataPdbAdj],
                 geoTriosA, title,printPath,fileName + "_A")
help.trioReports(["Unrestricted",dataPdbUn],
                 ["Restricted",dataPdbRes],
                 ["Reduced",dataPdbCut],
                 ["Adjusted",dataPdbAdj],
                 geoTriosB, title,printPath,fileName+ "_B")
'''
help.trioReports(["Unrestricted",dataPdbUn],
                 ["Restricted",dataPdbRes],
                 ["Reduced",dataPdbCut],
                 ["Adjusted",dataPdbAdj],
                 geoTriosTAU1, title,printPath,fileName+ "_t1")
geoTriosC = [
            ['PHI', 'PSI', 'C:O'],
            ['PHI', 'PSI', 'TAU'],
            ['PHI', 'PSI', 'TAU-1'],
            ['PHI', 'PSI', 'TAU+1'],
            ['PSI', 'N:N+1', 'C:O'],
            ['PSI', 'N:N+1', 'TAU'],
            ['PSI', 'N:N+1', 'TAU+1'],
            ['PSI', 'N:N+1', 'TAU-1'],
            ['TAU+1', 'TAU', 'C:O'],
            ['TAU','TAU-1', 'C:O'],
            ['TAU-1','TAU+1', 'C:O'],
           ]

'''
help.trioHexbins(["Unrestricted",dataPdbUn],
                 ["Restricted",dataPdbRes],
                 ["Reduced",dataPdbCut],
                 ["Adjusted",dataPdbAdj],
                 geoTriosC, title,printPath,fileName + "_hex")
'''



'''
help.trioReports(["Unrestricted",onePdbUn],
                 ["Restricted",onePdbRes],
                 ["Reduced",onePdbCut],
                 ["Adjusted",onePdbAdj],
                 geoTriosB, title,printPath,fileName + "_4r2x")
'''
help.trioReports(["Unrestricted",onePdbUn],
                 ["Restricted",onePdbRes],
                 ["Reduced",onePdbCut],
                 ["Adjusted",onePdbAdj],
                 geoTriosTAU1, title,printPath,fileName + "_4r2x_t1")
