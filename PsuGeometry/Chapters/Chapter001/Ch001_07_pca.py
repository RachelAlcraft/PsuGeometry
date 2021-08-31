'''
This file creates the data we want for PCA analysis. It creates in a dataframe which will be loaded into R
'''
import pandas as pd
from PsuGeometry import GeoReport as psu
import Ch000_Functions as help

filesPDBRoot ='C:/Dev/Github/ProteinDataFiles/pdb_data/'
filesADJRoot ='C:/Dev/Github/ProteinDataFiles/pdb_out/Fov2_ADJ/' #adjusted on Fo at 3 degrees thevenaz
loadPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/'
printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'

geos = ['PSI','PHI','N:CA','CA:C','C:N+1','C:O','N:CA:C','CA:C:O','O:C:N+1','CA:C:N+1','C:N+1:CA+1','C-1:N:CA','C-1:O-1','O:{,N,HOH,}','O:{,HOH,}','O:{,N,}']

fileName = 'PCA_CO'

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


print("---- Making adjusted CSV file for R PCA analysis--------")
dataPdbAdj = help.makeCsv('ADJUSTED', pdbListIn, geos, badAtoms,False)
dataPdbAdj = help.applyRestrictions(dataPdbAdj)
dataPdbAdj.to_csv(loadPath + "pca_co_adjusted.csv",index=False)