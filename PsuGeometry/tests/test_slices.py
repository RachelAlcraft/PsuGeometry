# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
'''
This script generates reports based on pure 2d mathematical matrices
The data could be anythong
It happens to be electron density of tyrosine rings
'''
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slicesgeo/'

filenames = ['Slices/1ejg/1TYR-1ejg.ccp4_slice.csv','Slices/1ejg/2TYR-1ejg.ccp4_slice.csv']
filenamesdiffs = ['Slices/1ejg/1TYR-1ejg_diff.ccp4_slice.csv','Slices/1ejg/2TYR-1ejg_diff.ccp4_slice.csv']
georep = psu.GeoReport([],pdbDataPath,edDataPath,printPath,ed=True,dssp=True)

slices = []
slicesdiffs = []

for fn in range(0,len(filenames)):
    slice = georep.loadSlice(filenames[fn])
    slicediff = georep.loadSlice(filenamesdiffs[fn])
    slices.append(slice)
    slicesdiffs.append(slicediff)

    georep.addSlice(slice, palette='cubehelix_r')
    georep.addSlice(slice, palette='cubehelix_r', logged=True)
    georep.addSlice(slicediff, palette='seismic', logged=False,centre=True)



georep.addSlices(slices, palette='cubehelix_r', title='Average', logged=False, centre=False)
georep.addSlices(slices, palette='cubehelix_r', title='Average', logged=True, centre=False)
georep.addSlices(slicesdiffs, palette='seismic', title='Average', logged=False, centre=True)


georep.printToHtml('Slices from Density Flight', 3, 'ccp4_slices')