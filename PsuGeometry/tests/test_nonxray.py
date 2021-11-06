
from PsuGeometry import GeoReport as psu

#pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
#edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
#printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Levels/'
#pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/LeicippusTesting/PdbFiles/'
edDataPath = 'C:/Dev/Github/ProteinDataFiles/ccp4_data/'
printPath = 'C:/Dev/Github/ProteinDataFiles/LeicippusTesting/'

pdbList = ['density_1ejg']#,'7a6a']#,'5xsg','6j60','6uos','6kj2','6bzm','6m9j','6cf4','6axz'] # high res cryoem
geoName = 'all'

# Create the geoemtric data
geos = ['PHI','PSI','N:CA','CA:C','C:O','C:N+1','TAU+1','N:N+1','TAU']
hueList = ['dssp','aa','bfactor','rid','atomNo'] # note the hues are the sum od the atoms

#Create the geoplots
for pdb in pdbList:
    georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath,ed=False,dssp=False)
    data = georep.getGeoemtryCsv(geos, hueList)
    data.to_csv(printPath + 'cryoem_'+pdb+'.csv', index=False)

    georep.addHistogram(data=data,geoX='N:CA',title='N-CA',ghost=True)
    georep.addHistogram(data=data,geoX='CA:C',title='CA-C',ghost=True)
    georep.addHistogram(data=data,geoX='C:O',title='C-O',ghost=True)
    georep.addHistogram(data=data,geoX='C:N+1',title='C-O',ghost=True)

    georep.addScatter(data=data,geoX='PHI',geoY='PSI',title='',hue='aa',palette='gist_rainbow',ghost=True)
    georep.addScatter(data=data, geoX='PHI', geoY='PSI', title='', hue='atomNo', palette='gist_rainbow', ghost=True)
    georep.addScatter(data=data, geoX='C:O', geoY='C:N+1', title='', hue='aa', palette='gist_rainbow', ghost=True)
    georep.addScatter(data=data, geoX='C:O', geoY='C:N+1', title='', hue='atomNo', palette='gist_rainbow', ghost=True)

    georep.addScatter(data=data, geoX='N:CA', geoY='CA:C', title='', hue='aa', palette='gist_rainbow', ghost=True)
    georep.addScatter(data=data, geoX='N:CA', geoY='CA:C', title='', hue='atomNo', palette='gist_rainbow', ghost=True)
    georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', title='', hue='aa', palette='gist_rainbow', ghost=True)
    georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', title='', hue='atomNo', palette='gist_rainbow', ghost=True)

    georep.addScatter(data=data, geoX='C:O', geoY='TAU+1', title='', hue='aa', palette='gist_rainbow', ghost=True)
    georep.addScatter(data=data, geoX='C:O', geoY='TAU+1', title='', hue='atomNo', palette='gist_rainbow', ghost=True)
    georep.addScatter(data=data, geoX='PSI', geoY='TAU', title='', hue='aa', palette='gist_rainbow', ghost=True)
    georep.addScatter(data=data, geoX='PSI', geoY='TAU', title='', hue='atomNo', palette='gist_rainbow', ghost=True)

    georep.printToHtml('Cryo EM Structures',4,'nonxray_' + pdb)
