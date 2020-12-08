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
#pdbList1000 = pdbList1000[:10]

georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False)

geoList = ['TAU', 'PSI']
hueList = ['aa', 'rid', 'resolution']  # note the hues are the sum of the atoms
# Create the dataframe
data = georep.getGeoemtryCsv(geoList, hueList)

print('Creating report')

for aa in ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']:
    georep.addDifference(dataA=data,dataB=data,geoX='N:CA:C', geoY='N:CA:C:N+1', restrictionsA={'aa': aa},exclusionsB={'aa': aa})

# Print the report
georep.printToHtml('Bimodal Tau', 3, 'Results2_tau')

