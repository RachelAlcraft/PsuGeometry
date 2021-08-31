###############################################
# This is to produce angle reports per amino acid include CIS PRO
###############################################
import pandas as pd
from PsuGeometry import GeoReport as psu
###############################################
fileUnrestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataUnrestricted.csv'
fileRestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataRestricted.csv'
fileReduced = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csv'
fileAdjusted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataAdjusted.csv'

#### HELPER FUNCTION ###############################

def evidenceReports(pdbSet,  fourSetNames, dataA, dataB, dataC, dataD ,trios, title,perAA=True, tag=''):
    pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'

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
                    georep.addScatter(data=dataA, geoX=trio[0], geoY=trio[1], hue=trio[2], title=aa + ' ' + trio[0] + '|' + trio[1] + '|' + trio[2] + ' ' + fourSetNames[0], palette='jet', sort='NON')
                    georep.addScatter(data=dataB, geoX=trio[0], geoY=trio[1], hue=trio[2], title=aa + ' ' + trio[0] + '|' + trio[1] + '|' + trio[2] + ' ' + fourSetNames[1], palette='jet',  sort='NON')
                    georep.addScatter(data=dataC, geoX=trio[0], geoY=trio[1], hue=trio[2],  title=aa + ' ' + trio[0] + '|' + trio[1] + '|' + trio[2] + ' ' + fourSetNames[2], palette='jet',  sort='NON')
                    georep.addScatter(data=dataD, geoX=trio[0], geoY=trio[1], hue=trio[2], title=aa + ' ' + trio[0] + '|' + trio[1] + '|' + trio[2] + ' ' + fourSetNames[3], palette='jet',  sort='NON')
                else:
                    georep.addHistogram(data=dataCutA, geoX=trio[0],title=aa + ' ' + fourSetNames[0] + ' ' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutB, geoX=trio[0],title=aa + ' ' + fourSetNames[1] + ' ' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutC, geoX=trio[0],title=aa + ' ' + fourSetNames[2] + ' ' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutD, geoX=trio[0],title=aa + ' ' + fourSetNames[3] + ' ' + trio[0], hue='ID')
        else:
            if len(trio) == 3:
                georep.addScatter(data=dataA, geoX=trio[0], geoY=trio[1], hue=trio[2],title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' ' + fourSetNames[0], palette='jet', sort='NON')
                georep.addScatter(data=dataB, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' ' + fourSetNames[1], palette='jet', sort='NON')
                georep.addScatter(data=dataC, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' ' + fourSetNames[2], palette='jet', sort='NON')
                georep.addScatter(data=dataD, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' ' + fourSetNames[3], palette='jet', sort='NON')
            else:
                georep.addHistogram(data=dataA, geoX=trio[0], title=fourSetNames[0] + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataB, geoX=trio[0], title=fourSetNames[1] + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataC, geoX=trio[0], title=fourSetNames[2] + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataD, geoX=trio[0], title=fourSetNames[3] + ' ' + trio[0], hue='ID')

    georep.printToHtml(title, 4, pdbSet + tag)

## PROGRAM EXECUTION ##################################################################
dataCsvAdjusted = pd.read_csv(fileAdjusted)
dataCsvReduced = pd.read_csv(fileReduced)
dataCsvRestricted = pd.read_csv(fileRestricted)
dataCsvUnrestricted = pd.read_csv(fileUnrestricted)

#manpulate omega values for planarity view
dataCsvAdjustedOmegaA = dataCsvAdjusted.query('OMEGA <= -90')
dataCsvAdjustedOmegaB = dataCsvAdjusted.query('OMEGA >= 90')
dataCsvAdjustedOmegaC = dataCsvAdjusted.query('OMEGA < 90')
dataCsvAdjustedOmegaC = dataCsvAdjustedOmegaC.query('OMEGA > -90')

dataCsvReducedOmegaA = dataCsvReduced.query('OMEGA <= -90')
dataCsvReducedOmegaB = dataCsvReduced.query('OMEGA >= 90')
dataCsvReducedOmegaC = dataCsvReduced.query('OMEGA < 90')
dataCsvReducedOmegaC = dataCsvReducedOmegaC.query('OMEGA > -90')

dataCsvRestrictedOmegaA = dataCsvRestricted.query('OMEGA <= -90')
dataCsvRestrictedOmegaB = dataCsvRestricted.query('OMEGA >= 90')
dataCsvRestrictedOmegaC = dataCsvRestricted.query('OMEGA < 90')
dataCsvRestrictedOmegaC = dataCsvRestrictedOmegaC.query('OMEGA > -90')

dataCsvUnrestrictedOmegaA = dataCsvUnrestricted.query('OMEGA <= -90')
dataCsvUnrestrictedOmegaB = dataCsvUnrestricted.query('OMEGA >= 90')
dataCsvUnrestrictedOmegaC = dataCsvUnrestricted.query('OMEGA < 90')
dataCsvUnrestrictedOmegaC = dataCsvUnrestrictedOmegaC.query('OMEGA > -90')

fileMergedDiffAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csvmean.csv'
fileUnrestrictedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataUnrestricted.csvmean.csv'
fileRestrictedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataRestricted.csvmean.csv'
fileReducedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csvmean.csv'
fileAdjustedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataAdjusted.csvmean.csv'

dataCsvAdjustedAv = pd.read_csv(fileAdjustedAv)
dataCsvReducedAv = pd.read_csv(fileMergedDiffAv)
dataCsvRestrictedAv = pd.read_csv(fileRestrictedAv)
dataCsvUnrestrictedAv = pd.read_csv(fileUnrestrictedAv)
dataCsvMergedAv = pd.read_csv(fileMergedDiffAv)
dataCsvAdjustedAv['ID'] = dataCsvAdjustedAv['pdbCode']
dataCsvReducedAv['ID'] = dataCsvReducedAv['pdbCode']
dataCsvRestrictedAv['ID'] = dataCsvRestrictedAv['pdbCode']
dataCsvUnrestrictedAv['ID'] = dataCsvUnrestrictedAv['pdbCode']
dataCsvMergedAv['ID'] = dataCsvMergedAv['pdbCode']


#broken into chunks for matplotlib memory allocation
chunkA, chunkB, chunkC,chunkD = False,False,False,True

if chunkA:
    geoTrios = [['N:CA'],['CA:C'],['C:O'],['C:N+1']]
    evidenceReports('Fo_Adj',['Unrestricted','Restricted','Reduced','Adjusted'],dataCsvUnrestricted,dataCsvRestricted,dataCsvReduced,dataCsvAdjusted,geoTrios,'Evidential Geometry Set Compare', perAA=False,tag='_bonds')
    geoTrios = [['TAU']]
    evidenceReports('Fo_Adj',['Unrestricted','Restricted','Reduced','Adjusted'],dataCsvUnrestricted,dataCsvRestricted,dataCsvReduced,dataCsvAdjusted,geoTrios,'Evidential Geometry Set Compare', perAA=True,tag='_tau')
    geoTrios = [['TAU-1']]
    evidenceReports('Fo_Adj',['Unrestricted','Restricted','Reduced','Adjusted'],dataCsvUnrestricted,dataCsvRestricted,dataCsvReduced,dataCsvAdjusted,geoTrios,'Evidential Geometry Set Compare', perAA=True,tag='_pretau')
    geoTrios = [['TAU+1']]
    evidenceReports('Fo_Adj',['Unrestricted','Restricted','Reduced','Adjusted'],dataCsvUnrestricted,dataCsvRestricted,dataCsvReduced,dataCsvAdjusted,geoTrios,'Evidential Geometry Set Compare', perAA=True,tag='_plustau')
if chunkB:
    geoTrios = [['CA:C:O']]
    evidenceReports('Fo_Adj',['Unrestricted','Restricted','Reduced','Adjusted'],dataCsvUnrestricted,dataCsvRestricted,dataCsvReduced,dataCsvAdjusted,geoTrios,'Evidential Geometry Set Compare', perAA=True,tag='_caco')
    geoTrios = [['O:C:N+1']]
    evidenceReports('Fo_Adj',['Unrestricted','Restricted','Reduced','Adjusted'],dataCsvUnrestricted,dataCsvRestricted,dataCsvReduced,dataCsvAdjusted,geoTrios,'Evidential Geometry Set Compare', perAA=True,tag='_ocn1')
    geoTrios = [['CA:C:N+1']]
    evidenceReports('Fo_Adj',['Unrestricted','Restricted','Reduced','Adjusted'],dataCsvUnrestricted,dataCsvRestricted,dataCsvReduced,dataCsvAdjusted,geoTrios,'Evidential Geometry Set Compare', perAA=True,tag='_cacn1')
if chunkC:
    geoTrios = [['OMEGA']]
    evidenceReports('Fo_Adj',['Unrestricted','Restricted','Reduced','Adjusted'],dataCsvUnrestrictedOmegaA,dataCsvRestrictedOmegaA,dataCsvReducedOmegaA,dataCsvAdjustedOmegaA,geoTrios,'Evidential Geometry Set Compare', perAA=False,tag='_omegaA')
    geoTrios = [['OMEGA']]
    evidenceReports('Fo_Adj',['Unrestricted','Restricted','Reduced','Adjusted'],dataCsvUnrestrictedOmegaB,dataCsvRestrictedOmegaB,dataCsvReducedOmegaB,dataCsvAdjustedOmegaB,geoTrios,'Evidential Geometry Set Compare', perAA=False,tag='_omegaB')
    geoTrios = [['OMEGA']]
    evidenceReports('Fo_Adj',['Unrestricted','Restricted','Reduced','Adjusted'],dataCsvUnrestrictedOmegaC,dataCsvRestrictedOmegaC,dataCsvReducedOmegaC,dataCsvAdjustedOmegaC,geoTrios,'Evidential Geometry Set Compare', perAA=False,tag='_omegaC')
if chunkD:
    geoTrios = [['N:CA'],['CA:C'],['C:O'],['C:N+1']]
    evidenceReports('Fo_Adj',['Unrestricted_Mean','Restricted_Mean','Reduced_Mean','Adjusted_Mean'],dataCsvUnrestrictedAv,dataCsvRestrictedAv,dataCsvReducedAv,dataCsvAdjustedAv,geoTrios,'Evidential Geometry Set Compare Averages', perAA=False,tag='_bondsMean')






