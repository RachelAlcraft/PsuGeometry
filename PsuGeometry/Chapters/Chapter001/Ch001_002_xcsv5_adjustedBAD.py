'''
This script creates a file for hand chosen geos and saves it
This particular set of geos is centred around the C:O question
'''

import pandas as pd
import Ch000_Functions as help
import matplotlib
print(matplotlib.__version__)

#filesPDBRoot ='C:/Dev/Github/ProteinDataFiles/pdb_data/'
#filesADJRoot ='C:/Dev/Github/ProteinDataFiles/pdb_out/Fov2_ADJ/' #adjusted on Fo at 3 degrees thevenaz
#loadPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/'
#printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'

geos = ['TAU','TAU+1','TAU-1','CA:C:O','O:C:N+1','CA-1:CA:CA+1',
        'N:CA:O','CA:O:N+1','O-1:N:CA',
        'O-1:C-1','C-1:N','N:CA','CA:C','C:O','C:N+1','N+1:CA+1','CA+1:C+1','C+1:O+1',
        'PHI','PSI','OMEGA','CA-1:C-1:N:CA',
        'CA-1:CA','CA:CA+1','C-1:C','C:C+1','N-1:N','N:N+1',
        'CA-1:N','CA-1:O-1','O-1:N','C-1:CA','N:C','CA:O','CA:N+1','O:N+1','C:CA+1','N+1:C+1',
        'O-1:CA','N:O','O:CA+1','N+1:O+1','N-1:O-1']


title='Backbone Report'
fileName = 'backbone'



print('### CREATING csv files ###')
pdbListIn = help.getPDBList()

#we want to look at ALL adjusted without the bad list

print("---- Making adjusted--------")
dataPdbAdj = help.makeCsv('ADJUSTED', pdbListIn, geos, [],False)
dataPdbAdj = help.applyRestrictions(dataPdbAdj)
dataPdbAdj = help.embellishCsv(dataPdbAdj)

# embellish with dssp - the dssp file was created ages ago from the linux laptop
pdbdssp = pd.read_csv('C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/CsvGeos_BEST_Set0DSSPALL.csv')
pdbdata = pd.read_csv('../../PdbLists/Pdbs_100.csv')

#embellish with dssp, resolution and software
try:
    pdbdssp['rid'] = pdbdssp['rid'].astype(str)
    pdbdssp['PDB'] =pdbdssp['pdbCode']
    pdbdssp['DSSPID'] = pdbdssp['pdbCode'] + pdbdssp['chain'] + pdbdssp['rid']
    pdbdssp = pdbdssp[['DSSPID', 'dssp']]
    pdbdata = pdbdata[['PDB', 'SOFTWARE', 'RES']]
    pdbdssp = pdbdssp.set_index('PDB').join(pdbdata.set_index('PDB'))
    pdbdssp['SOFTWARE'] = pdbdssp['SOFTWARE'].str[:8]
    pdbdssp = pdbdssp.dropna()

except:
    print('empty csv')

allList = []
allList.append([dataPdbAdj, help.loadPath + "bb_adjusted_BAD.csv"])

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