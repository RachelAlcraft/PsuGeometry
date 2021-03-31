# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
from PsuGeometry import Categoriser as cluster
'''
EH stats report comparison
'''
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
loadPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/'
printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/Reports/'

BestFileName = 'Set4_BestWithGeosALL.csv'
dataBest = pd.read_csv(loadPath + BestFileName)

produceRequests = True
aas = ['GLY']
geos = ['X']
'''
'N:O-2','C:O-2','N:CA:C:O-2','N:CA:N+1:O-2'
'''

#TIMER
print('----------start report 14----------')
startx = time.time()

for geo in geos:
    for aa in aas:
        bestCut = dataBest.query("aa ==  '" + aa + "'")
        bestCut['CLUSTER'] = 'X'  # cluster.tauCategory(data['PSI'],data['PHI'],data['N:N+1'],data['N:O-2'],data['N:CA:C:O-2'])
        # Apply my cluster function
        bestCut['CLUSTER'] = bestCut.apply(lambda row: cluster.tauCategory(row['PSI'], row['PHI'], row['N:N+1'], row['C:O-2'],row['C:{O}'],row['N:O-2'],row['N:CA:C:O-2'],row['N:CA:N+1:O-2'],row['N:{O}'],row['N:CA:C:{O}'],row['N:CA:N+1:{O}']), axis=1)
        cats = bestCut['CLUSTER'].unique()
        print(cats)
        for cat in cats:
            bestCutCat = bestCut.query("CLUSTER ==  '" + cat + "'")
            if produceRequests:
                bestCutCat['REQUEST'] = '?' # this prepares for the request to the electron density
                bestCutCat['REQUEST'] = bestCutCat.apply(lambda row: cluster.clusterEdTauMaker(row['pdbCode'], row['rid'], row['chain'], row['aa']), axis=1)
                saveCut = bestCutCat[['REQUEST']]
                saveCut.to_csv(loadPath + "Category_" + cat + ".csv", index=False,sep=' ')
                #And make a report for the category
            georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)
            georep.addScatter(data=bestCutCat, geoX='PHI', geoY='PSI', hue='TAU_x', title='', palette='jet')
            georep.addScatter(data=bestCutCat, geoX='PHI', geoY='PSI', hue='dssp', title='', palette='Dark2_r',categorical=True)
            georep.addScatter(data=bestCutCat, geoX='PSI', geoY='N:N+1', hue='TAU_x', title='', palette='jet')
            georep.addScatter(data=bestCutCat, geoX='PSI', geoY='N:N+1', hue='dssp', title='', palette='Dark2_r',categorical=True)
            georep.addScatter(data=bestCutCat, geoX='N:O-2', geoY='N:CA:C:O-2', hue='TAU_x', title='', palette='jet')
            georep.addScatter(data=bestCutCat, geoX='N:O-2', geoY='N:CA:C:O-2', hue='dssp', title='', palette='Dark2_r', categorical=True)
            georep.addScatter(data=bestCutCat, geoX='N:O-2', geoY='N:CA:N+1:O-2', hue='TAU_x', title='', palette='jet')
            georep.addScatter(data=bestCutCat, geoX='N:O-2', geoY='N:CA:N+1:O-2', hue='dssp', title='', palette='Dark2_r', categorical=True)
            georep.addScatter(data=bestCutCat, geoX='N:{O}', geoY='N:CA:C:{O}', hue='TAU_x', title='', palette='jet')
            georep.addScatter(data=bestCutCat, geoX='N:{O}', geoY='N:CA:C:{O}', hue='dssp', title='', palette='Dark2_r', categorical=True)
            georep.addScatter(data=bestCutCat, geoX='N:{O}', geoY='N:CA:N+1:{O}', hue='TAU_x', title='', palette='jet')
            georep.addScatter(data=bestCutCat, geoX='N:{O}', geoY='N:CA:N+1:{O}', hue='dssp', title='', palette='Dark2_r', categorical=True)
            title = 'DSSP key:H:a-helix S:bend G:3-helix E:extended strand \nT:h-bond turn B:isolated b-bridge I:5-helix'
            georep.printToHtml('Category ' + cat + '\n' + title, 4, 'BS_CAT_' + cat)





print('----------Finished----------')
endx = time.time()
time_diff = endx - startx
timestring = str(int(time_diff / 60)) + "m " + str(int(time_diff % 60)) + "s"
print(timestring)
