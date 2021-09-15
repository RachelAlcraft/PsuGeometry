'''
This script creates a file for hand chosen geos and saves it
This particular set of geos is centred around the C:O question
'''

import pandas as pd
import Ch000_Functions as help
import matplotlib
print(matplotlib.__version__)

geos = ['TAU','TAU+1','TAU-1','CA:C:O','O:C:N+1','CA-1:CA:CA+1',
        'N:CA:O','CA:O:N+1','O-1:N:CA',
        'O-1:C-1','C-1:N','N:CA','CA:C','C:O','C:N+1','N+1:CA+1','CA+1:C+1','C+1:O+1',
        'PHI','PSI','OMEGA','CA-1:C-1:N:CA',
        'CA-1:CA','CA:CA+1','C-1:C','C:C+1','N-1:N','N:N+1',
        'CA-1:N','CA-1:O-1','O-1:N','C-1:CA','N:C','CA:O','CA:N+1','O:N+1','C:CA+1','N+1:C+1',
        'O-1:CA','N:O','O:CA+1','N+1:O+1','N-1:O-1']

print('### CREATING csv files ###')
pdbListIn = help.getPDBList()
print("---- Getting bad atom list--------")
badAtoms = help.getBadAtomsListFromFile()  # Get the bad atoms list we will use to reduce the list further

print("---- Making adjusted--------")
dataPdbAdj = help.makeCsv('ADJUSTEDLAP', pdbListIn, geos, badAtoms,False)
#dataPdbAdj = pd.read_csv(help.loadPath + "bblap_adjusted_a.csv")
dataPdbAdj.to_csv(help.loadPath + "bblap_adjusted_a.csv", index=False)

dataPdbAdj = help.embellishCsv(dataPdbAdj)

dataPdbAdj = help.applyRestrictions(dataPdbAdj,True,True,True,True,False)
dataPdbAdj.to_csv(help.loadPath + "bblap_adjusted_b.csv", index=False)

print("---- Save to",help.loadPath + "bblap_adjusted.csv",'-------')
dataPdbAdj.to_csv(help.loadPath + "bblap_adjusted.csv", index=False)

