# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
import _Helpers as help
'''
Compare the maxima for a synthetic structure to the real structure
'''

def maximaCompareFake(pdbSet,pdbList,tag,reduce):
    print('Plotting fake maxima differences',pdbSet)

    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_out/' + pdbSet
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataM/'

    realCsv, badRealCsv, fakeCsv,badFakeCsv = help.getMaximaDiffs(pdbSet,pdbList,True)

    if reduce:
        #fakeCsv1 = fakeCsv.query('BGridDistance < 1.7')
        #fakeCsv2 = fakeCsv1.query('BGridDistance < 0.95')
        #fakeCsv1 = fakeCsv1.query('BGridDistance > 1.05')
        #fakeCsv3 = fakeCsv1.query('BGridDistance < 1.35')
        #fakeCsv1 = fakeCsv1.query('BGridDistance > 1.45')
        #fakeCsv = pd.concat([fakeCsv1,fakeCsv2,fakeCsv3])
        fakeCsv = fakeCsv.query('Difference <= 0.05')



        #fakeCsv2 = fakeCsv1.query('BGridDistance < 1.4 or BGridDistance > 1.42')
        #fakeCsv3 = fakeCsv2.query('BGridDistance < 0.95 or BGridDistance > 1.05')



    print(fakeCsv)



    fakeCol = 'jet_r'#'GnBu'

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

    georep.addHistogram(data=fakeCsv, geoX='Difference', title='Fake PDB structures, Occ=1, BFact=2', count=True, hue='pdbCode',palette='LightSeaGreen')
    georep.addScatter(data=fakeCsv, geoX='Difference', geoY='pdbCode', hue='Width', categorical=False,   palette=fakeCol + '', sort='RANDOM', title='')
    georep.addScatter(data=fakeCsv, geoX='Difference', geoY='Width', hue='pdbCode', categorical=True,   palette='tab20', sort='RANDOM', title='')

    georep.addScatter(data=fakeCsv, geoX='Difference', geoY='GridDistance', hue='BGridDistance', categorical=False, palette=fakeCol,  sort='RANDOM', title='')
    georep.addScatter(data=fakeCsv, geoX='Difference', geoY='BGridDistance', hue='GridDistance', categorical=False, palette=fakeCol, sort='RANDOM', title='')
    georep.addScatter(data=fakeCsv, geoX='GridDistance', geoY='BGridDistance', hue='Difference', categorical=False,  palette=fakeCol, sort='RANDOM', title='')

    georep.addProbability(data=fakeCsv, geoX='Difference', geoY='GridDistance',palette='cubehelix_r')
    georep.addProbability(data=fakeCsv, geoX='Difference', geoY='BGridDistance',palette='cubehelix_r')
    georep.addProbability(data=fakeCsv, geoX='GridDistance', geoY='BGridDistance',palette='cubehelix_r')


    for pdb in pdbList:
        print(pdb)
        georep.addHistogram(data=fakeCsv, geoX='Difference', title='Fake Density, Occ=1, BFact=10', hue='AtomNo', count=True, restrictions={'pdbCode': pdb}, palette='LightSeaGreen')
        georep.addHistogram(data=fakeCsv, geoX='Difference', title='Fake Density, Occ=1, BFact=10', hue='Reason', count=True, restrictions={'pdbCode': pdb}, palette='LightSeaGreen')
        georep.addScatter(data=fakeCsv, geoX='AtomNo', geoY='Difference', hue='AtomType', categorical=True,  restrictions={'pdbCode': pdb}, palette='tab20', sort='RANDOM', title='')



    georep.printToHtml('Maxima differences, set=' + pdbSet, 3, 'Maxima_' + pdbSet + tag)


def maximaCompareReal(pdbSet, pdbList, tag):
    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_out/' + pdbSet
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataM/'

    realCsv, badRealCsv, occRealCsv = help.getMaximaDiffs(pdbSet, pdbList, False)

    realOccOne = realCsv.query("Occupancy == 1")
    realCutDown = realOccOne.query("BFactor < 10")
    realCutDown5 = realCutDown.query("Difference <= 0.05")
    realCol = 'brg'#'RdPu'


    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

    georep.addHistogram(data=realCsv, geoX='Difference', title='PDB structures', count=True, hue='pdbCode')
    georep.addHistogram(data=realCutDown, geoX='Difference', title='PDB structures, Occ=1, BFact<10', count=True, hue='pdbCode')
    georep.addHistogram(data=realCutDown5, geoX='Difference', title='PDB structures, Occ=1, BFact<10, Difference<=0.05', count=True,   hue='pdbCode')
    georep.addScatter(data=realCsv, geoX='Difference', geoY='BFactor', hue='BGridDistance', palette=realCol, sort='RANDOM', title='Occ=1')
    georep.addScatter(data=realCutDown, geoX='Difference', geoY='BFactor', hue='BGridDistance', palette=realCol, sort='RANDOM', title='bfactor<=10')
    georep.addScatter(data=realCutDown5, geoX='Difference', geoY='BFactor', hue='BGridDistance', palette=realCol, sort='RANDOM', title='Diff<=0.05')
    georep.addScatter(data=realCsv, geoX='Difference', geoY='Width', hue='BFactor', categorical=False, palette=realCol, sort='RANDOM', title='Occ=1')
    georep.addScatter(data=realCutDown, geoX='Difference', geoY='Width', hue='BFactor', categorical=False, palette=realCol, sort='RANDOM', title='BFactor<=10')
    georep.addScatter(data=realCutDown5, geoX='Difference', geoY='Width', hue='BFactor', categorical=False, palette=realCol, sort='RANDOM', title='Diif<=0.05')

    #for pdb in pdbList:
    #    print(pdb)
    #    georep.addHistogram(data=realCsv, geoX='Difference', title='Occ=1', hue='AtomNo', count=True, restrictions={'pdbCode': pdb})
    #    georep.addHistogram(data=realCutDown, geoX='Difference', title='Occ=1, BFact<10', hue='AtomNo', count=True, restrictions={'pdbCode': pdb})
    #    georep.addScatter(data=realCutDown, geoX='Difference', geoY='BGridDistance', hue='GridDistance', categorical=False, restrictions={'pdbCode': pdb}, palette=realCol, sort='RANDOM', title='')
    #    georep.addScatter(data=realCutDown, geoX='Difference', geoY='AtomType', hue='AtomNo', categorical=False, restrictions={'pdbCode': pdb}, palette=realCol, sort='RANDOM', title='')


    georep.printToHtml('Maxima differences in PDB Structures, set=' + pdbSet, 3, 'Maxima_' + pdbSet + tag)



##########
# RUN TO TEST
#maximaCompare('CUBIC',help.getList('TOP20',0))

#maximaCompare('LINEAR',['3x2m','6e6o','4zm7','6s2m','1r6j','5yce','1x6z','5kwm','5d8v','2wfj','5nw3','1ucs','4rek','2b97'],'A')

#maximaCompare('LINEAR2',['3nir','5d8v','2wfj','5nw3','1ucs','3x2m','6e6o','4zm7','6s2m','1r6j','4rek','2b97','2ov0','2wfi','3zr8','5gji','3x34','5yce','1x6z','5kwm'],'B')
#maximaCompare('CUBIC2',['3nir','5d8v','2wfj','5nw3','1ucs','3x2m','6e6o','4zm7','6s2m','1r6j','4rek','2b97','2ov0','2wfi','3zr8','5gji','3x34','5yce','1x6z','5kwm'],'B')
#maximaCompare('QUINTIC2',['3nir','5d8v','2wfj','5nw3','1ucs','3x2m','6e6o','4zm7','6s2m','1r6j','4rek','2b97','2ov0','2wfi','3zr8','5gji','3x34','5yce','1x6z','5kwm'],'B')
#maximaCompare('HEPTIC2',['3nir','5d8v','2wfj','5nw3','1ucs','3x2m','6e6o','4zm7','6s2m','1r6j','4rek','2b97','2ov0','2wfi','3zr8','5gji','3x34','5yce','1x6z','5kwm'],'B')
#maximaCompare('NONIC2',['3nir','5d8v','2wfj','5nw3','1ucs','3x2m','6e6o','4zm7','6s2m','1r6j','4rek','2b97','2ov0','2wfi','3zr8','5gji','3x34','5yce','1x6z','5kwm'],'B')
#maximaCompare('HENDIC2',['3nir','5d8v','2wfj','5nw3','1ucs','3x2m','6e6o','4zm7','6s2m','1r6j','4rek','2b97','2ov0','2wfi'],'B') #,'3zr8','5gji','3x34','5yce','1x6z','5kwm'],'B')

'''
2wjf is a problem structure to investigate, no it was the wrong one!!!
'''
#maximaCompare('LINEAR',['2wjf'],'x')


def compareAtomsPdbAdjusted(dataCombined, geos, pdbSet,tag):
    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_out/' + pdbSet
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataM/'

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

    for geo in geos:

        georep.addHistogram(data=dataCombined, geoX=geo + '_Orig', title='Pdb Atoms', count=True, hue='pdbCode')
        georep.addHistogram(data=dataCombined, geoX=geo + '_Adj', title='Adjusted Atoms', count=True, hue='pdbCode')
        georep.addScatter(data=dataCombined, geoX=geo + '_Orig', geoY=geo + '_Adj', hue='SOFTWARE', palette='jet_r',sort='RANDOM', categorical=True,title='Comparing atom positions ' + geo)
        georep.addHexBins(data=dataCombined, geoX=geo + '_Orig', geoY=geo + '_Adj', title='Count ' + geo, hue='count', palette='cubehelix_r')
        #georep.addScatter(data=dataCombined, geoX=geo + '_Orig', geoY='SOFTWARE', hue=geo + '_Adj', palette='jet_r', sort='RANDOM', categorical=True, title='Comparing atom positions ' + geo)

    #for pdb in pdbList:
    #    print(pdb)
    #    georep.addHistogram(data=realCsv, geoX='Difference', title='Occ=1', hue='AtomNo', count=True, restrictions={'pdbCode': pdb})
    #    georep.addHistogram(data=realCutDown, geoX='Difference', title='Occ=1, BFact<10', hue='AtomNo', count=True, restrictions={'pdbCode': pdb})
    #    georep.addScatter(data=realCutDown, geoX='Difference', geoY='BGridDistance', hue='GridDistance', categorical=False, restrictions={'pdbCode': pdb}, palette=realCol, sort='RANDOM', title='')
    #    georep.addScatter(data=realCutDown, geoX='Difference', geoY='AtomType', hue='AtomNo', categorical=False, restrictions={'pdbCode': pdb}, palette=realCol, sort='RANDOM', title='')


    georep.printToHtml('Comparing atom positions: PDB vs maxima, set=' + pdbSet, 4, 'Compare_' + pdbSet + tag)



