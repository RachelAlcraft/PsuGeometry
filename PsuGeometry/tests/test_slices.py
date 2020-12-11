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

dir = 'Slices/1us0/'
fileroot = '1us0.ccp4_slice.csv'
filerootdiff = '1us0_diff.ccp4_slice.csv'
residues = ['A39','A48','A82','A102','A107','A177','A189','A198','A209','A291','A309']
pdb = '1us0'



georep = psu.GeoReport([],pdbDataPath,edDataPath,printPath,ed=True,dssp=True)

slices = []
slicesdiffs = []

for res in residues:
    slice = georep.loadSlice(dir+res+fileroot)
    slicediff = georep.loadSlice(dir+res+filerootdiff)
    slices.append(slice)
    slicesdiffs.append(slicediff)

    georep.addSlice(slice, palette='cubehelix_r',title=res)
    georep.addSlice(slice, palette='cubehelix_r', logged=True, title=res + " logged")
    georep.addSlice(slicediff, palette='seismic', logged=False,centre=True, title= res + ' difference')



georep.addSlices(slices, palette='cubehelix_r', title='Average', logged=False, centre=False)
georep.addSlices(slices, palette='cubehelix_r', title='Average logged', logged=True, centre=False)
georep.addSlices(slicesdiffs, palette='seismic', title='Average difference', logged=False, centre=True)


georep.printToHtml('Slices from Density Flight TYR for 1US0', 3, 'ccp4_slices' + pdb)