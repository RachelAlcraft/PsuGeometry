# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
import _Helpers as help
'''
EH stats report coparison
'''

def evidenceReports(pdbSet,  fourSetNames, dataA, dataB, dataC, dataD ,trios, title,perAA=True, tag=''):
    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_data/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataI/'

    aas = dataA['aa'].values
    aas = list(set(aas))
    aas.sort()

    georep = psu.GeoReport([],pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

    for trio in trios:
        if perAA:
            for aa in aas:
                dataCutA = dataA.query("aa ==  '" + aa + "'")
                dataCutB = dataB.query("aa ==  '" + aa + "'")
                dataCutC = dataC.query("aa ==  '" + aa + "'")
                dataCutD = dataD.query("aa ==  '" + aa + "'")
                if len(trio) == 3:
                    georep.addScatter(data=dataCutA, geoX=trio[0], geoY=trio[1], hue=trio[2], title=aa + ':' + trio[0] + ':'+ trio[1] , palette='jet', sort='NON')
                    georep.addScatter(data=dataCutB, geoX=trio[0], geoY=trio[1], hue=trio[2], title=aa + ':' + trio[0] + ':' + trio[1], palette='jet', sort='NON')
                    georep.addScatter(data=dataCutC, geoX=trio[0], geoY=trio[1], hue=trio[2],title=aa + ':' + trio[0] + ':' + trio[1], palette='jet', sort='NON')
                    georep.addScatter(data=dataCutD, geoX=trio[0], geoY=trio[1], hue=trio[2], title=aa + ':' + trio[0] + ':' + trio[1], palette='jet', sort='NON')
                else:
                    georep.addHistogram(data=dataCutA, geoX=trio[0],title=fourSetNames[0] + ' ' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutB, geoX=trio[0],title=fourSetNames[1] + ' ' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutC, geoX=trio[0],title=fourSetNames[2] + ' ' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutD, geoX=trio[0],title=fourSetNames[3] + ' ' + trio[0], hue='ID')
        else:
            if len(trio) == 3:
                georep.addScatter(data=dataA, geoX=trio[0], geoY=trio[1], hue=trio[2],title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' Unrestricted', palette='jet', sort='NON')
                georep.addScatter(data=dataB, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' Restricted', palette='jet', sort='NON')
                georep.addScatter(data=dataC, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' Restricted+cut', palette='jet', sort='NON')
                georep.addScatter(data=dataD, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' Adjusted', palette='jet', sort='NON')
            else:
                georep.addHistogram(data=dataA, geoX=trio[0], title=fourSetNames[0] + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataB, geoX=trio[0], title=fourSetNames[1] + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataC, geoX=trio[0], title=fourSetNames[2] + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataD, geoX=trio[0], title=fourSetNames[3] + ' ' + trio[0], hue='ID')



    georep.printToHtml(title, 4, pdbSet + '_Defensible' + tag)
