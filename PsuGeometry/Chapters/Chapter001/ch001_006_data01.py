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

title='Hydrogen bonding Report'
fileName = 'hydrogenbonding'


print('### LOADING csv files ###') # bit rubbish but we didn;t change the object references with dssp
#dataPdbUn = pd.read_csv(help.loadPath + "hb_unrestricted.csv")
#dataPdbRes = pd.read_csv(help.loadPath + "hb_restricted.csv")
dataPdbCut = pd.read_csv(help.loadPath + "hb_reduced.csv")
dataPdbAdj = pd.read_csv(help.loadPath + "hb_adjusted.csv")

# ensure data is correctly restricted
#dataPdbUn = help.applyRestrictions(dataPdbUn,True,False,False,False)
dataPdbCut = help.applyRestrictions(dataPdbCut,True,True,True,True)
#dataPdbRes = help.applyRestrictions(dataPdbRes,True,True,True,True)
dataPdbAdj = help.applyRestrictions(dataPdbAdj,True,True,True,False)

tag = ''
#SHale we cut on bfactor factor?
BFactorFactor = False
if BFactorFactor:
    tag = '_bff'
    dataPdbCut = dataPdbCut.query('bfactorRatio <= 1.2')
    #dataPdbRes = dataPdbRes.query('bfactorRatio <= 1.2')
    dataPdbAdj = dataPdbAdj.query('bfactorRatio <= 1.2')


print('### Creating scatter files ###')
'''
geos = ['TAU','TAU+1',
        'N:CA','CA:C','C:O','C:N+1',
        'CA:C:O:N+1',
        'O:{,HOH,}','N:{,HOH,}',
        'O:{,O,}','N:{,O,}',
        'O:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}','N:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}']
'''

geoTriosA = [
            ['C:O','O:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}','NearO',True],
            ['C:O','N+1:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}','NearN+1',True],
            ['C:O','O:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}','RES',False],
            ['C:O','N+1:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}','RES',False],
            ['C:O','O:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}','SOFTWARE',True],
            ['C:O','N+1:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}','SOFTWARE',True],
            ['O:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}','N+1:{,O,}','C:O',False],
            ['C:O','LapDiff', 'RES',False],
            ['C:O','O:{,HOH,}', 'bfactor',False],
            ['C:O','O:{,O,}', 'bfactor',False],
            ['C:O','O:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}', 'bfactor',False],
            ['C:O','N+1:{,HOH,}', 'bfactor',False],
            ['C:O','N+1:{,O,}', 'bfactor',False],
            ['C:O','N+1:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}', 'bfactor',False],
           ]

help.trioReports(
                 ["Reduced",dataPdbCut],
                 ["Adjusted",dataPdbAdj],
                 ["NONE",""],
                 ["NONE",""],
                 geoTriosA, title,help.printPath,fileName + tag)

