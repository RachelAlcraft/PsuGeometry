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
dataPdbCut = pd.read_csv(help.loadPath + "bb_reduced01.csv")
dataPdbAdj = pd.read_csv(help.loadPath + "bb_adjusted01.csv")

# Prepare a seperate report for the pdb
onePdbUn = dataPdbUn.query("pdbCode == '4r2x'")
onePdbRes = dataPdbRes.query("pdbCode == '4r2x'")
onePdbCut = dataPdbCut.query("pdbCode == '4r2x'")
onePdbAdj = dataPdbAdj.query("pdbCode == '4r2x'")

#Choose a bfactor to cut at
#dataPdbUn = dataPdbUn.query("`C:O_avbfactor` <= 15")
#dataPdbRes = dataPdbRes.query("`C:O_avbfactor` <= 15")
#dataPdbCut = dataPdbCut.query("`C:O_avbfactor` <= 15")
#dataPdbAdj = dataPdbAdj.query("`C:O_avbfactor` <= 15")

dataPdbUn = dataPdbUn.query("pdbCode != '4r2x'")
dataPdbRes = dataPdbRes.query("pdbCode != '4r2x'")
dataPdbCut = dataPdbCut.query("pdbCode != '4r2x'")
dataPdbAdj = dataPdbAdj.query("pdbCode != '4r2x'")


print('### Creating scatter files ###')

geoTriosA = [['C:O'],
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
            ['TAU+1', 'TAU', 'C:O',False],
            ['TAU','TAU-1', 'C:O',False],
            ['TAU-1','TAU+1', 'C:O',False],
            ['C:O','TAU+1', 'bfactor',False],
           ]

help.trioReports(["Unrestricted",dataPdbUn],
                 ["Restricted",dataPdbRes],
                 ["Reduced",dataPdbCut],
                 ["Adjusted",dataPdbAdj],
                 geoTriosA, title,help.printPath,fileName + "_a01")

help.trioReports(["Unrestricted",onePdbUn],
                 ["Restricted",onePdbRes],
                 ["Reduced",onePdbCut],
                 ["Adjusted",onePdbAdj],
                 geoTriosA, title,help.printPath,fileName + "_a01_4r2x")
