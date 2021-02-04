# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
'''
Proof of bimodal tau
'''
myWindowsLaptop = True
###############################################################################################
pdbList1000 = geol.GeoPdbLists().getListPaper()
pdbList1000 = pdbList1000[:100]

geos = ['TAU','PSI','N:O','N:N+1','O-1:N+1','CA-1:CA:CA+1:CA+2','CA-1:N+1','N:CA']
aas = ['GLY','ALA']
#plots = ['hex','scatter','count','probability','dssp']
plots = ['scatter']

dsspHue = 'dssp'
hueList = ['aa', 'rid', 'bfactor','dssp']

countPalette = 'cubehelix_r'
scatterPalette = 'jet'
hexPalette = 'jet'
dsspPalette = 'gist_ncar'
probPalette = 'ocean_r'

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

georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=False)

# Create the dataframe
dataX = georep.getGeoemtryCsv(geos, hueList)

geos = geos[1:]

#data1 = dataX.query('TAU <= 100')
#data2 = dataX.query('TAU > 100')

datas = [dataX]
tags = ['']
for i in range(0,len(datas)):
    data = datas[i]
    tag = tags[i]
    for aa in aas:
        sql = 'aa == "' + aa + '"'
        dataAll = data.query(sql)

        usedList = []

        for geoA in geos:
            for geoB in geos:
                if geoA != geoB:
                    if geoA + geoB not in usedList and geoB+geoA not in usedList:
                        for plot in plots:
                            if plot == 'hex':
                                georep.addHexBins(data=dataAll, geoX=geoA, geoY=geoB, hue='TAU', title=geoA + ':' + geoB,palette=hexPalette,bins=bins, gridsize=gridsize)
                            elif plot == 'count':
                                georep.addHexBins(data=dataAll, geoX=geoA, geoY=geoB, hue='count', title=geoA + ':' + geoB, palette=countPalette,bins=bins, gridsize=gridsize)
                            elif plot == 'scatter':
                                georep.addScatter(data=dataAll, geoX=geoA, geoY=geoB, hue='TAU', title=geoA + ':' + geoB,palette=scatterPalette, sort='NON')
                            elif plot == 'probability':
                                georep.addProbability(geoX=geoA, geoY=geoB, title=geoA + ':' + geoB, palette=probPalette)
                            elif plot == 'dssp':
                                georep.addScatter(data=dataAll, geoX=geoA, geoY=geoB, hue=dsspHue, title=geoA + ':' + geoB,palette=dsspPalette, sort='NON')

                            usedList.append(geoA + geoB)

        print('Creating reports')
        georep.printToHtml('Multi-modal Tau ' + aa + ' Plots ' + tag, 4, 'Results2_' + aa + '_taucorr_' + tag)





