
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/1000Structures/'

from PsuGeometry import GeoReport as geor
import pandas as pd

pdbList = []
pdbListPath = printPath + 'pro Structures.csv'
pdbdata = pd.read_csv(pdbListPath)
pdbList = pdbdata['pdb_code']
pdbList = pdbList[0:50]

georep = geor.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=False)

geoList = ['PHI','PSI','TAU','C-1:C','C-1:N:CA','CHI1','CA-1:CA','OMEGA','CA:CA+1','C-1:N:CA:C']
hueList = ['aa','bfactor','rid','resolution']

data = georep.getGeoemtryCsv(geoList, hueList)



georep.addScatter(data=data, geoX='PHI',geoY='PSI',hue='CHI1', title='Ramachandran', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='PHI',geoY='PSI',hue='CA-1:CA', title='Ramachandran', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})

georep.addScatter(data=data, geoX='OMEGA',geoY='PHI',hue='CHI1', title='Omega-Phi', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='OMEGA',geoY='PHI',hue='CA-1:CA', title='Omega-Phi', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})

georep.addScatter(data=data, geoX='C-1:C',geoY='PHI',hue='CHI1', title='C-1:C/Phi', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='C-1:C',geoY='PHI',hue='CA-1:CA', title='C-1:C/Phi', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})

georep.addScatter(data=data, geoX='C-1:N:CA',geoY='TAU',hue='CHI1', title='Tau angles', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='C-1:N:CA',geoY='TAU',hue='CA-1:CA', title='Tau angles', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})

georep.addScatter(data=data, geoX='C-1:C',geoY='C-1:N:CA:C',hue='CHI1', title='C-1:C/Phi', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='C-1:C',geoY='C-1:N:CA:C',hue='CA-1:CA', title='C-1:C/Phi', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})



# And finally create the reort with a file name of choice
georep.printToHtml('Correlations with geometric hue',2,'geohue')
