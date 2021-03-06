# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
'''
Proof of bimodal tau
'''
myWindowsLaptop = False
###############################################################################################
pdbList1000 = geol.GeoPdbLists().getListPaper()
pdbList1000 = pdbList1000[:500]
#pdbList1000 = ['5zj8']

dihs = ['PSI','PHI']
angles = ['TAU']
distances = ['N:O'] # OCA


aas = ['GLY','ALA']

dsspHue = 'dssp'
hueList = ['aa', 'rid', 'bfactor']

palette1 = 'viridis'
palette2 = 'rainbow'
bins=50
gridsize=50

###################################################################################
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'
includeDSSP = False
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper01/'
    includeDSSP = False  # on my windows computer
if not includeDSSP:
    hueList = ['aa', 'rid', 'bfactor']
    dsspHue = 'aa'
###########################################################################################





georepData = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False)

geoList = []
for geo in dihs:
    geoList.append(geo)
for geo in distances:
    geoList.append(geo)
for geo in angles:
    geoList.append(geo)

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

        #try:
        #    georep.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='PSI', restrictionsA={'aa': aa},exclusionsB={'aa': aa})
        #    georep.addHistogram(geoX='TAU', data=dataAll, title=aa, hue=dsspHue)
        #    georep.addHexBins(data=dataAll, geoX='TAU', geoY='PSI', hue='PHI', palette=palette2, bins=bins,gridsize=gridsize)
        #except:
        #    print('No data')

        georepData.addHexBins(data=dataAll, geoX='PHI', geoY='PSI', hue='count', palette=palette1, bins=bins, gridsize=gridsize)
        georepData.addHexBins(data=dataAll, geoX='PHI', geoY='PSI', hue='TAU', palette=palette2, bins=bins,gridsize=gridsize)
        georepData.addScatter(data=dataLower, geoX='PHI', geoY='PSI', hue='TAU', title='Rama Lower', palette=palette2, sort='NON')
        georepData.addScatter(data=dataMiddle, geoX='PHI', geoY='PSI', hue='TAU', title='Rama Middle', palette=palette2, sort='NON')
        georepData.addScatter(data=dataUpper, geoX='PHI', geoY='PSI', hue='TAU', title='Rama Upper', palette=palette2, sort='NON')

        hus = ['PSI','TAU']
        xaxs = ['TAU','PSI']
        pcs = ['plasma','magma']
        pss = ['jet', 'rainbow']

        for geo in geoList:
            for i in range(0,2):
                hu = hus[i]
                xax = xaxs[i]
                p1 = pcs[i]
                p2 = pss[i]

                georepData.addHexBins(data=dataAll, geoX=xax, geoY=geo, hue='count', palette=p1,bins=bins, gridsize=gridsize)
                if geo == 'PSI' and hu=='PSI':
                    hu = 'PHI'
                elif geo == 'TAU' and hu=='TAU':
                    hu = 'PHI'
                georepData.addHexBins(data=dataAll, geoX=xax, geoY=geo, hue=hu, palette=p2, bins=bins,gridsize=gridsize)
                georepData.addScatter(data=dataLower, geoX=xax, geoY=geo, hue=hu, title=geo + ' Lower', palette=p2, sort='NON')
                georepData.addScatter(data=dataMiddle, geoX=xax, geoY=geo, hue=hu, title=geo + ' Middle', palette=p2, sort='NON')
                georepData.addScatter(data=dataUpper, geoX=xax, geoY=geo, hue=hu, title=geo + ' Upper', palette=p2, sort='NON')

        print('Creating reports')
        georepData.printToHtml('Tau Correlations ' + aa + ' Plots, Pdbs=' + tag + str(len(pdbList1000)) , 5, 'Results2_' + aa + '_tauscratch_' + tag + str(len(pdbList1000)))





