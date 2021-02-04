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
#pdbList1000 = pdbList1000[:100]

dihs = ['PSI','PHI','CA-2:CA-1:CA:CA+1','CA-1:CA:CA+1:CA+2']
distances = ['CA-1:N+1','N:CA','CA:C','N:C','O-1:O','CA-1:CA','CA:CA+1','N-1:N','N:N+1','N:O','O-1:N','O:O+1','O-1:O+1','O-1:N+1']
angles = ['N:CA:C','CA:C:O']
aas = ['GLY','ALA']

hueList = ['aa', 'rid', 'bfactor','dssp']

palette1 = 'cubehelix_r'
palette2 = 'jet'
palette3 = 'gist_ncar'
bins='log'
gridsize=50

###################################################################################






georepData = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)

geoList = []
for geo in dihs:
    geoList.append((geo))
for geo in distances:
    geoList.append((geo))
for geo in angles:
    geoList.append((geo))

# Create the dataframe
dataAll = georepData.getGeoemtryCsv(geoList, hueList)

data1 = dataAll.query('TAU <= 100')
data2 = dataAll.query('TAU > 100')
data3 = dataAll.query('PSI<-100')
data4 = dataAll.query('PSI>100')
data5 = dataAll.query('PSI>=-100 and PSI <=100')

datas = [data2,data3,data4,data5]
tags = ['all','lower','upper','middle']
for i in range(0,len(datas)):
    data = datas[i]
    tag = tags[i]
    for aa in aas:
        sql = 'aa == "' + aa + '"'
        dataa = data.query(sql)
        georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)

        try:
            georep.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='PSI', restrictionsA={'aa': aa},exclusionsB={'aa': aa})
            georep.addHistogram(geoX='TAU', data=dataa, title=aa, hue='dssp')
        except:
            print('No data')

        georep.addHexBins(data=dataa, geoX='PHI', geoY='PSI', hue='count', palette=palette1,bins=bins, gridsize=gridsize)
        georep.addHexBins(data=dataa, geoX='PHI', geoY='PSI', hue='TAU', palette=palette2,bins=bins, gridsize=gridsize)
        georep.addScatter(data=dataa, geoX='PHI', geoY='PSI', hue='TAU', title=geo, palette=palette2, sort='NON')
        georep.addScatter(data=dataa, geoX='PHI', geoY='PSI', hue='dssp', title=geo, palette=palette3, sort='NON')


        for geo in dihs:
            georep.addHexBins(data=dataa, geoX='TAU', geoY=geo, hue='count', palette=palette1,bins=bins, gridsize=gridsize)
            if geo == 'PSI':
                georep.addHexBins(data=dataa, geoX='TAU', geoY=geo, hue='PHI', palette=palette2,bins=bins, gridsize=gridsize)
                georep.addScatter(data=dataa, geoX='TAU', geoY=geo, hue='PHI', title=geo, palette=palette2, sort='NON')
            else:
                georep.addHexBins(data=dataa, geoX='TAU', geoY=geo, hue='PSI', palette=palette2,bins=bins, gridsize=gridsize)
                georep.addScatter(data=dataa, geoX='TAU', geoY=geo, hue='PSI', title=geo, palette=palette2, sort='NON')
            georep.addScatter(data=dataa, geoX='TAU', geoY=geo, hue='dssp', title=geo, palette=palette3, sort='NON')


        for geo in angles:
            georep.addHexBins(data=dataa, geoX='TAU', geoY=geo, hue='count', palette=palette1,bins=bins, gridsize=gridsize)
            georep.addHexBins(data=dataa, geoX='TAU', geoY=geo, hue='PSI', palette=palette2,bins=bins, gridsize=gridsize)
            georep.addScatter(data=dataa, geoX='TAU', geoY=geo, hue='PSI', title=geo, palette=palette2, sort='NON')
            georep.addScatter(data=dataa, geoX='TAU', geoY=geo, hue='dssp', title=geo, palette=palette3, sort='NON')


        for geo in distances:
            georep.addHexBins(data=dataa, geoX='TAU', geoY=geo, hue='count', palette=palette1,bins=bins, gridsize=gridsize)
            georep.addHexBins(data=dataa, geoX='TAU', geoY=geo, hue='PSI', palette=palette2,bins=bins, gridsize=gridsize)
            georep.addScatter(data=dataa, geoX='TAU', geoY=geo, hue='PSI', title=geo, palette=palette2, sort='NON')
            georep.addScatter(data=dataa, geoX='TAU', geoY=geo, hue='dssp', title=geo, palette=palette3, sort='NON')

        print('Creating reports')
        georep.printToHtml('Multi-modal Tau ' + aa + ' Plots ' + tag, 4, 'Results2_' + aa + '_tau_' + tag)





