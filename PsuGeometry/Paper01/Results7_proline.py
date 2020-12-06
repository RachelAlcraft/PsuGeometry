# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
'''
This script looks at the distribution structural feautures of proline
By looking at some plots with geoemtric hues
'''


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

pdbList1000 = geol.GeoPdbLists().getList1000()
geoList = ['PHI','PSI','TAU','C-1:C','C-1:N:CA','CHI1','CA-1:CA','OMEGA','CA:CA+1','C-1:N:CA:C']
hueList = ['aa','bfactor','rid','resolution','pdbCode']

georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False)
data = georep.getGeoemtryCsv(geoList,hueList)

data = data[data['C-1:C'] < 10]


georep.addScatter(data=data, geoX='PHI',geoY='PSI',hue='CHI1', title='Ramachandran', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='PHI',geoY='PSI',hue='CA-1:CA', title='Ramachandran', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='OMEGA',geoY='PHI',hue='CHI1', title='Omega-Phi', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='OMEGA',geoY='PHI',hue='CA-1:CA', title='Omega-Phi', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='C-1:C',geoY='PHI',hue='CHI1', title='C-1:C/Phi', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='C-1:C',geoY='PHI',hue='CA-1:CA', title='C-1:C/Phi', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='C-1:N:CA',geoY='TAU',hue='CHI1', title='Tau angles', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='C-1:N:CA',geoY='TAU',hue='CA-1:CA', title='Tau angles', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})

title = 'Correlations with geometric hue<br/> Proline up-pucker by proxy of CHI1<br/>Proline cis/trans by proxy of short Calpha'
georep.printToHtml(title,2,'Results7_proline')