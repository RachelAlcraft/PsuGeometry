
'''
In this file we compare individual geos
to see if any pdbs are problematic
'''

import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoReport as psu

print('### LOADING csv files ###')
dataPdbUn = pd.read_csv(help.loadPath + "bb_unrestricted.csv")
dataPdbRes = pd.read_csv(help.loadPath + "bb_restricted.csv")
dataPdbCut = pd.read_csv(help.loadPath + "bb_reduced.csv")
dataPdbAdj = pd.read_csv(help.loadPath + "bbden_adjusted.csv")
dataPdbLap = pd.read_csv(help.loadPath + "bblap_adjusted.csv")
# ensure data is correctly restricted
dataPdbUn = help.applyRestrictions(dataPdbUn,True,False,False,False,False)
dataPdbRes = help.applyRestrictions(dataPdbRes,True,True,True,True,False)
dataPdbCut = help.applyRestrictions(dataPdbCut,True,True,True,True,True)
dataPdbAdj = help.applyRestrictions(dataPdbAdj,True,True,True,False,True)
dataPdbLap = help.applyRestrictions(dataPdbLap,True,True,True,False,True)

tag = ''
#SHale we cut on bfactor factor?
BFactorFactor = True
if BFactorFactor:
    tag = '_bff'
    dataPdbRes = dataPdbRes.query('bfactorRatio <= 1.2')
    dataPdbCut = dataPdbCut.query('bfactorRatio <= 1.2')
    dataPdbAdj = dataPdbAdj.query('bfactorRatio <= 1.2')

dsspList = dataPdbUn["dssp"].unique()

georepCO = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)
georepN1 = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

georepCO_comp = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)
georepN1_comp = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)


print('### DSSP reports ###')
georepCO.addHistogram(data=dataPdbUn, geoX='C:O',title='Unrestricted all', hue='ID')
georepCO.addHistogram(data=dataPdbRes, geoX='C:O',title='Restricted all', hue='ID')
georepCO.addHistogram(data=dataPdbCut, geoX='C:O',title='Cut 05 all', hue='ID')
georepCO.addHistogram(data=dataPdbAdj, geoX='C:O',title='Adjusted 05 all', hue='ID')

georepN1.addHistogram(data=dataPdbUn, geoX='C:N+1',title='Unrestricted all', hue='ID')
georepN1.addHistogram(data=dataPdbRes, geoX='C:N+1',title='Restricted all', hue='ID')
georepN1.addHistogram(data=dataPdbCut, geoX='C:N+1',title='Cut 05 all', hue='ID')
georepN1.addHistogram(data=dataPdbAdj, geoX='C:N+1',title='Adjusted 05 all', hue='ID')

for dssp in dsspList:
    print('### ',dssp,' ###')
    try:
        onepdbUn = dataPdbUn.query("dssp == '" + dssp + "'")
        onepdbRes = dataPdbRes.query("dssp == '" + dssp + "'")
        onepdbCut = dataPdbCut.query("dssp == '" + dssp + "'")
        onepdbAdj = dataPdbAdj.query("dssp == '" + dssp + "'")

        georepCO.addHistogram(data=onepdbUn, geoX='C:O', title='Unrestricted ' + str(dssp), hue='ID')
        georepCO.addHistogram(data=onepdbRes, geoX='C:O', title='Restricted ' + str(dssp), hue='ID')
        georepCO.addHistogram(data=onepdbCut, geoX='C:O', title='Cut 05 ' + str(dssp), hue='ID')
        georepCO.addHistogram(data=onepdbAdj, geoX='C:O', title='Adjusted 05 ' + str(dssp), hue='ID')

        georepN1.addHistogram(data=onepdbUn, geoX='C:N+1', title='Unrestricted ' + str(dssp), hue='ID')
        georepN1.addHistogram(data=onepdbRes, geoX='C:N+1', title='Restricted ' + str(dssp), hue='ID')
        georepN1.addHistogram(data=onepdbCut, geoX='C:N+1', title='Cut 05 ' + str(dssp), hue='ID')
        georepN1.addHistogram(data=onepdbAdj, geoX='C:N+1', title='Adjusted 05 ' + str(dssp), hue='ID')

        georepCO_comp.addHistogram(data=onepdbAdj, geoX='C:O', title='Adjusted ' + str(dssp), hue='ID')
        georepCO_comp.addStatsCompare(dataA=onepdbAdj, dataB=dataPdbAdj,descA=dssp,descB="All",geoX="C:O")
        georepCO_comp.addHistogram(data=dataPdbAdj, geoX='C:O', title='Adjusted All', hue='ID')

        georepN1_comp.addHistogram(data=onepdbAdj, geoX='C:N+1', title='Adjusted ' + str(dssp), hue='ID')
        georepN1_comp.addStatsCompare(dataA=onepdbAdj, dataB=dataPdbAdj, descA=dssp, descB="All", geoX="C:N+1")
        georepN1_comp.addHistogram(data=dataPdbAdj, geoX='C:N+1', title='Adjusted All', hue='ID')




    except:
        print('Error with the dssp code',dssp)

georepCO.printToHtml('DSSP Per C:O', 4, 'dssp_co' + tag)
georepN1.printToHtml('DSSP Per C:N+1', 4, 'dssp_n1'+tag)
georepCO_comp.printToHtml('DSSP against All, C:O', 3, 'dssp_co_comp'+tag)
georepN1_comp.printToHtml('DSSP against All, C:N+1', 3, 'dssp_n1_comp'+tag)



