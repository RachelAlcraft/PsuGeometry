

'''
This script creates a file for hand chosen geos and saves it
This particular set of geos is centred around the C:O question
'''

import pandas as pd
import Ch000_Functions as help
import matplotlib
print(matplotlib.__version__)

geos = ['TAU','TAU+1',
        'N:CA','CA:C','C:O','C:N+1',
        'CA:C:O:N+1',
        'O:{,HOH,}','N+1:{,HOH,}',
        'O:{,O,}','N+1:{,O,}',
        'O:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}','N+1:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}',
        'O:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}','N+1:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}']

print('### CREATING cut csv file ###')
pdbListIn = help.getPDBList()
#pdbListIn =pdbListIn[:10]
print("---- Getting bad atom list--------")
badAtoms = help.getBadAtomsListFromFile()  # Get the bad atoms list we will use to reduce the list further

print("---- Making reduced--------")
dataPdbCut = help.makeCsv('PDB', pdbListIn, geos, badAtoms,False)
#dataPdbCut = pd.read_csv(help.loadPath + "hb_reduced_a.csv")
dataPdbCut.to_csv(help.loadPath + "hb_reduced_a.csv", index=False)
dataPdbCut = help.applyRestrictions(dataPdbCut,True,True,True,True)
dataPdbCut.to_csv(help.loadPath + "hb_reduced_b.csv", index=False)
dataPdbCut = help.embellishCsv(dataPdbCut)

#make new columns for hues
dataPdbCut['NearN+1'] = dataPdbCut['N+1:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}_atmmotif']
dataPdbCut['NearO'] = dataPdbCut['O:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}_atmmotif']

print("---- Save to",help.loadPath + "hb_reduced.csv",'-------')
dataPdbCut.to_csv(help.loadPath + "hb_reduced.csv", index=False)

print("---- Making adjusted--------")
dataPdbAdj = help.makeCsv('ADJUSTED', pdbListIn, geos, badAtoms,False)
#dataPdbAdj = pd.read_csv(help.loadPath + "hb_adjusted_a.csv")
dataPdbAdj.to_csv(help.loadPath + "hb_adjusted_a.csv", index=False)
dataPdbAdj = help.applyRestrictions(dataPdbAdj,True,True,True,True)
dataPdbAdj.to_csv(help.loadPath + "hb_adjusted_b.csv", index=False)
dataPdbAdj = help.embellishCsv(dataPdbAdj)

#make new columns for hues
dataPdbAdj['NearN+1'] = dataPdbAdj['N+1:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}_atmmotif']
dataPdbAdj['NearO'] = dataPdbAdj['O:{,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,O,HOH,}_atmmotif']

print("---- Save to",help.loadPath + "hb_adjusted.csv",'-------')
dataPdbAdj.to_csv(help.loadPath + "hb_adjusted.csv", index=False)

