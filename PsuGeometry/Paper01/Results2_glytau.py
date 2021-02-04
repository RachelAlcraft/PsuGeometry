# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
'''
Proof of bimodal tau
'''
myWindowsLaptop = True
###############################################################################################
pdbList1000 = geol.GeoPdbLists().getListPaper()
pdbList1000 = pdbList1000[50:100]

dihs = ['PSI','PHI','OMEGA','CA-2:CA-1:CA:CA+1','CA-1:CA:CA+1:CA+2','N:CA:C:O','O:C:N+1:CA+1']
distances = ['CA-1:N+1','N:CA','CA:C','N:C','O-1:O','CA-1:CA','CA:CA+1','N-1:N','N:N+1','N:O','O-1:N','O:O+1','O-1:O+1','O-1:N+1','C-1:N','C:N+1','O-1:N+2']
angles = ['N:CA:C','CA:C:O','C-1:N:CA','CA:C:N+1','CA-1:CA:C','N-1:N:N+1','C-1:C:C+1','O-1:O:O+1']
aas = ['GLY','ALA']

dsspHue = 'dssp'
hueList = ['aa', 'rid', 'bfactor','dssp']

palette1 = 'cubehelix_r'
palette2 = 'jet'
palette3 = 'gist_ncar'
bins='log'
gridsize=50

includeDSSP = True
###################################################################################
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper01/'
    includeDSSP = False  # on my windows computer
if not includeDSSP:
    hueList = ['aa', 'rid', 'bfactor']
    dsspHue = 'aa'
###########################################################################################





georepData = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=False)

geoList = []
for geo in dihs:
    geoList.append((geo))
for geo in distances:
    geoList.append((geo))
for geo in angles:
    geoList.append((geo))

# Create the dataframe
dataX = georepData.getGeoemtryCsv(geoList, hueList)

data1 = dataX.query('TAU <= 100')
data2 = dataX.query('TAU > 100')

datas = [data2]
tags = ['']
for i in range(0,len(datas)):
    data = datas[i]
    tag = tags[i]
    for aa in aas:
        sql = 'aa == "' + aa + '"'
        dataAll = data.query(sql)

        dataLower = dataAll.query('PSI<-100')
        dataMiddle = dataAll.query('PSI>=-100 and PSI <=100')
        dataUpper = dataAll.query('PSI>100')

        georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)

        try:
            georep.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='PSI', restrictionsA={'aa': aa},exclusionsB={'aa': aa})
            georep.addHistogram(geoX='TAU', data=dataAll, title=aa, hue=dsspHue)
        except:
            print('No data')

        georep.addHexBins(data=dataAll, geoX='PHI', geoY='PSI', hue='count', palette=palette1,bins=bins, gridsize=gridsize)
        georep.addHexBins(data=dataLower, geoX='PHI', geoY='PSI', hue='count', palette=palette1, bins=bins,gridsize=gridsize)
        georep.addHexBins(data=dataMiddle, geoX='PHI', geoY='PSI', hue='count', palette=palette1, bins=bins,gridsize=gridsize)
        georep.addHexBins(data=dataUpper, geoX='PHI', geoY='PSI', hue='count', palette=palette1, bins=bins,gridsize=gridsize)

        for geo in dihs:
            georep.addHexBins(data=dataAll, geoX='TAU', geoY=geo, hue='count', palette=palette1,bins=bins, gridsize=gridsize)
            hu = 'PSI'
            if geo == 'PSI':
                hu = 'PHI'
            georep.addScatter(data=dataLower, geoX='TAU', geoY=geo, hue=hu, title=geo + ' Lower', palette=palette2, sort='NON')
            georep.addScatter(data=dataMiddle, geoX='TAU', geoY=geo, hue=hu, title=geo + ' Middle', palette=palette2, sort='NON')
            georep.addScatter(data=dataUpper, geoX='TAU', geoY=geo, hue=hu, title=geo + ' Upper', palette=palette2, sort='NON')


        for geo in angles:
            georep.addHexBins(data=dataAll, geoX='TAU', geoY=geo, hue='count', palette=palette1,bins=bins, gridsize=gridsize)
            georep.addScatter(data=dataLower, geoX='TAU', geoY=geo, hue='PSI', title=geo + ' Lower', palette=palette2,sort='NON')
            georep.addScatter(data=dataMiddle, geoX='TAU', geoY=geo, hue='PSI', title=geo + ' Middle', palette=palette2,sort='NON')
            georep.addScatter(data=dataUpper, geoX='TAU', geoY=geo, hue='PSI', title=geo + ' Upper', palette=palette2,sort='NON')


        for geo in distances:
            georep.addHexBins(data=dataAll, geoX='TAU', geoY=geo, hue='count', palette=palette1, bins=bins,gridsize=gridsize)
            georep.addScatter(data=dataLower, geoX='TAU', geoY=geo, hue='PSI', title=geo + ' Lower', palette=palette2,sort='NON')
            georep.addScatter(data=dataMiddle, geoX='TAU', geoY=geo, hue='PSI', title=geo + ' Middle', palette=palette2,sort='NON')
            georep.addScatter(data=dataUpper, geoX='TAU', geoY=geo, hue='PSI', title=geo + ' Upper', palette=palette2,sort='NON')

        print('Creating reports')
        georep.printToHtml('Multi-modal Tau ' + aa + ' Plots ' + tag, 4, 'Results2_' + aa + '_tau_' + tag)





