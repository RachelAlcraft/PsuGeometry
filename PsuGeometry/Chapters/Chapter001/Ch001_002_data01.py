'''
This script creates a file for hand chosen geos, so it is slower as they are no serialised.
This particular report is designed to look at hydrogen bonding on the carbonyl oxygen
'''

import pandas as pd
from PsuGeometry import GeoReport as psu
import Ch000_Functions as help
import matplotlib
print(matplotlib.__version__)

#filesPDBRoot ='C:/Dev/Github/ProteinDataFiles/pdb_data/'
#filesADJRoot ='C:/Dev/Github/ProteinDataFiles/pdb_out/Fov2_ADJ/' #adjusted on Fo at 3 degrees thevenaz
#loadPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/'
#printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'

geos = ['TAU','TAU+1','TAU-1','CA:C:O','O:C:N+1','CA:C:N+1','CA-1:CA:CA+1',
        'O-1:C-1','C-1:N','N:CA','CA:C','C:O','C:N+1','N+1:CA+1','CA+1:C+1','C+1:O+1',
        'PHI','PSI','OMEGA','CA-1:C-1:N:CA',
        'CA-1:CA','CA:CA+1','C-1:C','C:C+1','N-1:N','N:N+1',
        'CA-1:N','CA-1:O-1','O-1:N','C-1:CA','N:C','CA:O','CA:N+1','O:N+1','C:CA+1','N+1:C+1',
        'O-1:CA','N:O','O:CA+1','N+1:O+1','N-1:O-1']


title='Backbone Report'
fileName = 'backbone'


print('### LOADING csv files ###') # bit rubbish but we didn;t change the object references with dssp
dataPdbUn = pd.read_csv(help.loadPath + "bb_unrestricted.csv")
dataPdbRes = pd.read_csv(help.loadPath + "bb_restricted.csv")
dataPdbCut = pd.read_csv(help.loadPath + "bb_reduced.csv")
dataPdbAdj = pd.read_csv(help.loadPath + "bb_adjusted.csv")

# ensure data is correctly restricted
dataPdbUn = help.applyRestrictions(dataPdbUn,True,False,False,False)
dataPdbCut = help.applyRestrictions(dataPdbCut,True,True,True,True)
dataPdbRes = help.applyRestrictions(dataPdbRes,True,True,True,True)
dataPdbAdj = help.applyRestrictions(dataPdbAdj,True,True,True,False)

tag = '_a'
#SHale we cut on bfactor factor?
BFactorFactor = True
if BFactorFactor:
    tag = '_bff'
    dataPdbCut = dataPdbCut.query('bfactorRatio <= 1.2')
    dataPdbRes = dataPdbRes.query('bfactorRatio <= 1.2')
    dataPdbAdj = dataPdbAdj.query('bfactorRatio <= 1.2')


print('### Creating scatter files ###')

geoTriosA = [['N:CA'],['CA:C'],['C:O'],['C:N+1'],
            ['CA:O:N+1', 'TAU+1', 'C:O', False],
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
            ['PSI', 'C:N+1', 'C:O',False],
            ['PSI', 'C:O', 'C:N+1',False],
            ['C:O', 'C:N+1', 'PSI',False],
            ['C:O', 'C:N+1', 'OMEGA',False],
            ['C:O', 'C:N+1', 'CA-1:C-1:N:CA',False],
            ['TAU+1', 'TAU', 'C:O',False],
            ['TAU','TAU-1', 'C:O',False],
            ['TAU-1','TAU+1', 'C:O',False],
            ['C:O','TAU+1', 'bfactor',False],
            ['C:O','C:N+1', 'TAU',False],
            ['C:O','C:N+1', 'TAU+1',False],
            ['C:O','C:N+1', 'bfactor',False],
            ['C:O','C:N+1', 'dssp',True],
            ['C:O','bfactor', 'TAU+1',False],
            ['C:O','bfactor', 'COUNT',False],
           ]

help.trioReports(["Unrestricted",dataPdbUn],
                 ["Restricted",dataPdbRes],
                 ["Reduced",dataPdbCut],
                 ["Adjusted",dataPdbAdj],
                 geoTriosA, title,help.printPath,fileName + tag)

