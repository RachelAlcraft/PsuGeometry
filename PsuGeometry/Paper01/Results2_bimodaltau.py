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

A = False # Big scatter 5
B = False # Histograms
C = False # Angle scatter
D = False # Phi prob
E = True # psi prob
F = False # omega prob
G = False # nc prob
H = False # Phi/psi/tau rotation



georepA = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
georepB = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
georepC = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
georepD = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
georepE = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
georepF = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
georepG = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
georepH = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
#georepI = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)


geoList = ['TAU', 'PSI','PHI','OMEGA','N:O','N:C','N:CA','CA:C','CA:C:O','CA:O','C-1:N:CA','CA-1:CA','CA:C:N+1']
hueList = ['aa', 'rid', 'bfactor','dssp']  # note the hues are the sum of the atoms
# Create the dataframe
data = georepA.getGeoemtryCsv(geoList, hueList)

if A:
    georepA.addScatter(data=data, geoX='N:CA:C', geoY='aa', hue='dssp', title='Tau Data', palette='gist_ncar', sort='NON')
    georepA.addScatter(data=data, geoX='N:CA:C', geoY='PSI', hue='dssp', title='Tau/Psi', palette='gist_ncar', sort='NON')
    georepA.addScatter(data=data, geoX='PSI', geoY='aa', hue='dssp', title='Psi Data', palette='gist_ncar', sort='NON')
    georepA.addScatter(data=data, geoX='N:CA:C', geoY='N:O', hue='dssp', title='Tau/N:O', palette='gist_ncar', sort='NON')
    georepA.addScatter(data=data, geoX='N:CA:C', geoY='N:C', hue='dssp', title='Tau/N:C', palette='gist_ncar', sort='NON')
if H:
    georepH.addScatter(data=data, geoX='PHI', geoY='PSI', hue='TAU', title='PHI/PSI/TAU', palette='jet', sort='NON',vmin=95,vmax=135)
    georepH.addScatter(data=data, geoX='TAU', geoY='PSI', hue='PHI', title='PHI/PSI/TAU', palette='jet', sort='NON',vmin=-180,vmax=180)
    georepH.addScatter(data=data, geoX='TAU', geoY='PHI', hue='PSI', title='PHI/PSI/TAU', palette='jet', sort='NON',vmin=-180,vmax=180)
    georepH.addScatter(data=data, geoX='TAU', geoY='PHI', hue='bfactor', title='PHI/PSI/TAU', palette='cubehelix_r', sort='ASC')
    georepH.addScatter(data=data, geoX='PHI', geoY='PSI', hue='TAU', title='PHI/PSI/TAU', palette='jet', sort='NON',exclusions={'aa': 'GLY,PRO'},vmin=95,vmax=135)
    georepH.addScatter(data=data, geoX='TAU', geoY='PSI', hue='PHI', title='PHI/PSI/TAU', palette='jet', sort='NON',exclusions={'aa': 'GLY,PRO'},vmin=-180,vmax=180)
    georepH.addScatter(data=data, geoX='TAU', geoY='PHI', hue='PSI', title='PHI/PSI/TAU', palette='jet', sort='NON',exclusions={'aa': 'GLY,PRO'},vmin=-180,vmax=180)
    georepH.addScatter(data=data, geoX='TAU', geoY='PHI', hue='bfactor', title='PHI/PSI/TAU', palette='cubehelix_r', sort='ASC',exclusions={'aa': 'GLY,PRO'})
if C:
    georepC.addScatter(data=data, geoX='TAU', geoY='C-1:N:CA', hue='CA-1:CA', title='Backbone angles-1', palette='jet', sort='NON',vmin=2.8,vmax=4.2)
    georepC.addScatter(data=data, geoX='TAU', geoY='C-1:N:CA', hue='PHI', title='Backbone angles-1', palette='jet',sort='NON',vmin=-180,vmax=180)
    georepC.addScatter(data=data, geoX='TAU', geoY='C-1:N:CA', hue='PSI', title='Backbone angles-1', palette='jet', sort='NON',vmin=-180,vmax=180)
    georepC.addScatter(data=data, geoX='TAU', geoY='CA:C:N+1', hue='PHI', title='Backbone angles+1', palette='jet',sort='NON', vmin=-180, vmax=180)
    georepC.addScatter(data=data, geoX='TAU', geoY='CA:C:N+1', hue='PSI', title='Backbone angles+1', palette='jet',sort='NON', vmin=-180, vmax=180)

#for aa in ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']:
for aa in ['GLY']:
    print(aa)
    if A:
        georepA.addScatter(data=data, geoX='N:CA:C', geoY='PHI', hue='dssp', title='Phi/Tau',restrictions={'aa': aa}, palette='gist_ncar',sort='NON')
        georepA.addScatter(data=data, geoX='N:CA:C', geoY='PSI', hue='dssp', title='Psi/Tau',restrictions={'aa': aa}, palette='gist_ncar',sort='NON')
        georepA.addScatter(data=data, geoX='N:CA:C', geoY='OMEGA', hue='dssp', title='Omega/Tau', restrictions={'aa': aa}, palette='gist_ncar',sort='NON')
        georepA.addScatter(data=data, geoX='N:C', geoY='N:O', hue='dssp', title='N:C/N:O', restrictions={'aa': aa}, palette='gist_ncar',sort='NON')
        georepA.addScatter(data=data, geoX='N:CA:C', geoY='N:C', hue='dssp', title='Tau/N:C', restrictions={'aa': aa}, palette='gist_ncar',sort='NON')
    if C:
        georepC.addScatter(data=data, geoX='TAU', geoY='C-1:N:CA', hue='CA-1:CA', title='Backbone angles-1', palette='jet', sort='NON',restrictions={'aa': aa},vmin=2.8,vmax=4.2)
        georepC.addScatter(data=data, geoX='TAU', geoY='C-1:N:CA', hue='PHI', title='Backbone angles-1', palette='jet', sort='NON', restrictions={'aa': aa},vmin=-180,vmax=180)
        georepC.addScatter(data=data, geoX='TAU', geoY='C-1:N:CA', hue='PSI', title='Backbone angles-1', palette='jet', sort='NON', restrictions={'aa': aa},vmin=-180,vmax=180)
        georepC.addScatter(data=data, geoX='TAU', geoY='CA:C:N+1', hue='PHI', title='Backbone angles+!', palette='jet', sort='NON', restrictions={'aa': aa}, vmin=-180, vmax=180)
        georepC.addScatter(data=data, geoX='TAU', geoY='CA:C:N+1', hue='PSI', title='Backbone angles+!', palette='jet', sort='NON', restrictions={'aa': aa}, vmin=-180, vmax=180)

    if H:
        georepH.addScatter(data=data, geoX='PHI', geoY='PSI', hue='TAU', title='PHI/PSI/TAU', restrictions={'aa': aa}, palette='jet', sort='NON',vmin=95,vmax=135)
        georepH.addScatter(data=data, geoX='TAU', geoY='PSI', hue='PHI', title='PHI/PSI/TAU', restrictions={'aa': aa}, palette='jet', sort='NON',vmin=-180,vmax=180)
        georepH.addScatter(data=data, geoX='TAU', geoY='PHI', hue='PSI', title='PHI/PSI/TAU', restrictions={'aa': aa}, palette='jet', sort='NON',vmin=-180,vmax=180)
        georepH.addScatter(data=data, geoX='TAU', geoY='PHI', hue='bfactor', title='PHI/PSI/TAU', palette='cubehelix_r',sort='ASC', restrictions={'aa': aa})
    if B:
        georepB.addHistogram(geoX='N:CA:C',data=data,title=aa,restrictions={'aa': aa},hue='dssp')

    exc = aa
    if aa not in ('PRO,GLY'):
        exc += ',PRO,GLY'

    if D:
        georepD.addDifference(dataA=data, dataB=data, geoX='N:CA:C', geoY='PHI', restrictionsA={'aa': aa},exclusionsB={'aa': exc})
    if E:
        georepE.addDifference(dataA=data, dataB=data, geoX='N:CA:C', geoY='PSI', restrictionsA={'aa': aa},exclusionsB={'aa': exc})
    if F:
        georepF.addDifference(dataA=data, dataB=data, geoX='N:CA:C', geoY='OMEGA', restrictionsA={'aa': aa},exclusionsB={'aa': exc})
    if G:
        georepG.addDifference(dataA=data, dataB=data, geoX='N:CA:C', geoY='N:C', restrictionsA={'aa': aa},exclusionsB={'aa': exc})




print('Creating reports')

# Print the report
if A:
    georepA.printToHtml('Multi-modal Tau Scatter Plots', 5, 'Results2_tau_scatter')
if B:
    georepB.printToHtml('Tau Histograms', 4, 'Results2_tau_histograms')
if C:
    georepC.printToHtml('Tau Backbone Angles', 5, 'Results2_tau_angles')
if D:
    georepD.printToHtml('Multi-modal Tau PHI', 3, 'Results2_tau_phi_prob')
if E:
    georepE.printToHtml('Multi-modal Tau PSI', 3, 'Results2_tau_psi_prob')
if F:
    georepF.printToHtml('Multi-modal Tau OMEGA', 3, 'Results2_tau_omega_prob')
if G:
    georepG.printToHtml('Multi-modal Tau N:C', 3, 'Results2_tau_nc_prob')
if H:
    georepH.printToHtml('Multi-modal Tau Geoemtric Hues', 4, 'Results2_tau_phi_psi')





