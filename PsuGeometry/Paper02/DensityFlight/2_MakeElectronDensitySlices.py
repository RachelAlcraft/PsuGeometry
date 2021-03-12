# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import random
import pandas as pd
'''
TAU correlations
'''
###############################################################################################
myWindowsLaptop = True
FileDir = 'Set3B'
#Tag = 'Extreme values about PSI=0'
#Tag = 'Wide N, PSI +ve near 180'
Tag = 'Wide N, PSI -ve near 180'
#Tag = 'Short N, PSI < 5'
#Tag = 'Short N, PSI > 10'

###################################################################################
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper02/DensityFlight/'
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    edSlicePath = 'F:/Code/ProteinDataFiles/ccp4_out/'
    #printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'



#We are going to load the data that has been created by density flight
edSlicePath += FileDir + "/"
print(edSlicePath + "_Results.csv")
inputdata = pd.read_csv(edSlicePath + "_Results.csv")
#PdbCode,Tag,
# CentreX,CentreY,CentreZ,# LinearX,LinearY,LinearZ,PlanarX,PlanarY,PlanarZ,
# CentreV,LinearV,PlanarV,Angle,
# BCentreX,BCentreY,BCentreZ,BLinearX,BLinearY,BLinearZ,BPlanarX,BPlanarY,BPlanarZ,
# BCentreC,BLinearV,BPlanarV,BAngle

georep = psu.GeoReport([],pdbDataPath,edDataPath,edSlicePath,ed=False,dssp=False)

pdbs = inputdata['PdbCode'].values
tags = inputdata['Tag'].values
taus = inputdata['Angle'].values
btaus = inputdata['BAngle'].values

origs = []
betters = []
radiants = []
brads = []


# Once the app has created the data we can load it
for i in range(0,len(pdbs)):
    pdb = pdbs[i]
    tag = tags[i]
    tau = taus[i]
    btau = btaus[i]

    sliceOrigVal = georep.loadSlice(edSlicePath + pdb + tag + "value_slice.csv")
    sliceBetterVal = georep.loadSlice(edSlicePath + pdb + tag + "bvalue_slice.csv")
    sliceOrigRad = georep.loadSlice(edSlicePath + pdb + tag + "radiant_slice.csv")
    sliceBetterRad = georep.loadSlice(edSlicePath + pdb + tag + "bradiant_slice.csv")

    georep.addSlice(sliceOrigVal, palette='cubehelix_r',title=pdb + tag + ' value, tau=' + str(round(tau,3)))
    georep.addSlice(sliceBetterVal, palette='cubehelix_r',title=pdb + tag + ' better value, tau=' + str(round(btau,3)))
    georep.addSlice(sliceOrigRad, palette='bone',title=pdb + tag + ' radiant',Contour=False)
    georep.addSlice(sliceBetterRad, palette='bone',title=pdb + tag + ' better radiant',Contour=False)

    origs.append(sliceOrigVal)
    betters.append(sliceBetterVal)
    radiants.append(sliceOrigRad)
    brads.append(sliceBetterRad)

georep.addSlices(origs, palette='cubehelix_r', title='Average values', logged=False, centre=False)
georep.addSlices(betters, palette='cubehelix_r', title='Average better values', logged=False, centre=False)
georep.addSlices(radiants, palette='bone', title='Average radiant', logged=False, centre=False,Contour=False)
georep.addSlices(brads, palette='bone', title='Average better radiant', logged=False, centre=False,Contour=False)


georep.printToHtml(Tag,4,'_Results_ccp4_slices')