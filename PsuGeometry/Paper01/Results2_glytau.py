# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
'''
Proof of bimodal tau
'''


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

###############################################################################################

pdbList1000 = geol.GeoPdbLists().getListPaper()
pdbList1000 = pdbList1000[:100]



georepGLY = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
georepALA = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)

geoList = ['TAU', 'PSI','PHI','OMEGA','N:O','N:C','N:CA','CA:C','CA:C:O','CA:O','C-1:N:CA','CA-1:CA','CA:C:N+1','CA-2:CA-1:CA:CA+1','N-1:N','N:N+1','CA-1:CA','CA:CA+1']
hueList = ['aa', 'rid', 'bfactor','dssp']  # note the hues are the sum of the atoms
# Create the dataframe
data = georepALA.getGeoemtryCsv(geoList, hueList)
datagly = data.query('aa =="GLY"')
dataala = data.query('aa =="ALA"')

georepALA.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='PSI', restrictionsA={'aa': 'ALA'},exclusionsB={'aa': 'ALA,GLY'})
georepGLY.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='PSI', restrictionsA={'aa': 'GLY'},exclusionsB={'aa': 'GLY'})

georepALA.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='CA:C:O', restrictionsA={'aa': 'ALA'},exclusionsB={'aa': 'ALA,GLY'})
georepGLY.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='CA:C:O', restrictionsA={'aa': 'GLY'},exclusionsB={'aa': 'GLY'})

georepALA.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='CA-1:CA', restrictionsA={'aa': 'ALA'},exclusionsB={'aa': 'ALA,GLY'})
georepGLY.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='CA-1:CA', restrictionsA={'aa': 'GLY'},exclusionsB={'aa': 'GLY'})

georepALA.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='CA:CA+1', restrictionsA={'aa': 'ALA'},exclusionsB={'aa': 'ALA,GLY'})
georepGLY.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='CA:CA+!', restrictionsA={'aa': 'GLY'},exclusionsB={'aa': 'GLY'})

georepALA.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='CA-2:CA-1:CA:CA+1', restrictionsA={'aa': 'ALA'},exclusionsB={'aa': 'ALA,GLY'})
georepGLY.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='CA-2:CA-1:CA:CA+1', restrictionsA={'aa': 'GLY'},exclusionsB={'aa': 'GLY'})

georepALA.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='N:O', restrictionsA={'aa': 'ALA'},exclusionsB={'aa': 'ALA,GLY'})
georepGLY.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='N:O', restrictionsA={'aa': 'GLY'},exclusionsB={'aa': 'GLY'})

georepALA.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='N-1:N', restrictionsA={'aa': 'ALA'},exclusionsB={'aa': 'ALA,GLY'})
georepGLY.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='N-1:N', restrictionsA={'aa': 'GLY'},exclusionsB={'aa': 'GLY'})

georepALA.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='N:N+1', restrictionsA={'aa': 'ALA'},exclusionsB={'aa': 'ALA,GLY'})
georepGLY.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='N:N+1', restrictionsA={'aa': 'GLY'},exclusionsB={'aa': 'GLY'})

georepALA.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='N:CA', restrictionsA={'aa': 'ALA'},exclusionsB={'aa': 'ALA,GLY'})
georepGLY.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='N:CA', restrictionsA={'aa': 'GLY'},exclusionsB={'aa': 'GLY'})

georepALA.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='CA:C', restrictionsA={'aa': 'ALA'},exclusionsB={'aa': 'ALA,GLY'})
georepGLY.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='CA:C', restrictionsA={'aa': 'GLY'},exclusionsB={'aa': 'GLY'})

georepALA.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='N:C', restrictionsA={'aa': 'ALA'},exclusionsB={'aa': 'ALA,GLY'})
georepGLY.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='N:C', restrictionsA={'aa': 'GLY'},exclusionsB={'aa': 'GLY'})


georepALA.addHistogram(geoX='N:CA:C',data=dataala,title='TAU')
georepGLY.addHistogram(geoX='N:CA:C',data=datagly,title='TAU')

georepALA.addScatter(data=dataala, geoX='N:CA:C', geoY='PHI', hue='PSI', title='Tau/Phi', palette='jet', sort='NON')
georepGLY.addScatter(data=datagly, geoX='N:CA:C', geoY='PHI', hue='PSI', title='Tau/Phi', palette='jet', sort='NON')

georepALA.addScatter(data=dataala, geoX='N:CA:C', geoY='PHI', hue='dssp', title='Tau/Phi', palette='gist_ncar', sort='NON')
georepGLY.addScatter(data=datagly, geoX='N:CA:C', geoY='PHI', hue='dssp', title='Tau/Phi', palette='gist_ncar', sort='NON')

print('Creating reports')
georepALA.printToHtml('Multi-modal Tau ALA Plots', 3, 'Results2_ala_tau')
georepGLY.printToHtml('Multi-modal Tau Gly Plots', 3, 'Results2_gly_tau')





