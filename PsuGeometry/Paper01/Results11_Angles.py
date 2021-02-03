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

A = True  # Scatter
B = True  # Histograms
C = True  # Cross Correlations

georepA = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
georepB = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
georepC = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)


geoListA = ['CB:CA:C','N:CA:CB','N:CA:C','CA:C:O','CA:C:N+1','O:C:N+1','C-1:N:CA','PHI','PSI']
geoListB = ['N:CA:C','CA:C:O','CA:C:N+1','O:C:N+1','C-1:N:CA','PHI','PSI'] # for glycine as no CB
geoListGly = ['N:CA:C','CA:C:O','CA:C:N+1','O:C:N+1','C-1:N:CA']
geoListNonGly = ['CB:CA:C','N:CA:CB'] # for glycine as no CB

hueList = ['aa', 'rid', 'bfactor', 'dssp']  # note the hues are the sum of the atoms
# Create the dataframe
dataA = georepA.getGeoemtryCsv(geoListA, hueList)
dataB = georepA.getGeoemtryCsv(geoListB, hueList)

if A:
    for geo in geoListNonGly:
        georepA.addScatter(data=dataA, geoX='PHI', geoY='PSI', hue=geo, title=geo, palette='jet',sort='NON')
    for geo in geoListGly:
        georepA.addScatter(data=dataB, geoX='PHI', geoY='PSI', hue=geo, title=geo, palette='jet',sort='NON')

if B:
    for geo in geoListNonGly:
        georepB.addHistogram(geoX=geo,data=dataA,title=geo)
    for geo in geoListGly:
        georepB.addHistogram(geoX=geo,data=dataB,title=geo)
if C:
    for geoa in geoListGly:
        for geob in geoListGly:
            if geoa != geob:
                georepC.addScatter(data=dataB, geoX=geoa, geoY=geob, hue='PSI', title=geoa + ':' + geob, palette='jet',sort='NON')
        for geob in geoListNonGly:
            georepC.addScatter(data=dataA, geoX=geoa, geoY=geob, hue='PSI', title=geoa + ':' + geob, palette='jet',sort='NON')


for aa in ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']:
    print(aa)
    if A:
        for geo in geoListNonGly:
            georepA.addScatter(data=dataA, geoX='PHI', geoY='PSI', hue=geo, title=aa + ':' + geo, palette='jet', sort='NON',restrictions={'aa': aa})
        for geo in geoListGly:
            georepA.addScatter(data=dataB, geoX='PHI', geoY='PSI', hue=geo, title=aa + ':' + geo, palette='jet', sort='NON',restrictions={'aa': aa})

    if B:
        for geo in geoListNonGly:
            georepB.addHistogram(geoX=geo,data=dataA,title=aa + ':' + geo,restrictions={'aa': aa})
        for geo in geoListGly:
            georepB.addHistogram(geoX=geo,data=dataB,title=aa + ':' + geo,restrictions={'aa': aa})

    if C:
        for geoa in geoListGly:
            for geob in geoListGly:
                if geoa != geob:
                    georepC.addScatter(data=dataB, geoX=geoa, geoY=geob, hue='PSI', title=aa + ':' + geoa + ':' + geob, palette='jet',sort='NON',restrictions={'aa': aa})
            for geob in geoListNonGly:
                georepC.addScatter(data=dataA, geoX=geoa, geoY=geob, hue='PSI', title=aa + ':' + geoa + ':' + geob, palette='jet',sort='NON',restrictions={'aa': aa})

print('Creating reports')

# Print the report
if A:
    georepA.printToHtml('Angle Ramachandrans', 7, 'Results11_angle_ramas')
if B:
    georepB.printToHtml('Angle Histograms', 7, 'Results11_angle_histograms')
if C:
    georepC.printToHtml('Angle Scatters', 4, 'Results11_angle_scatters')





