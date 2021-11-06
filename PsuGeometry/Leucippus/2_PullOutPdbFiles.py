import FilesAndCSVs as fac
from PsuGeometry import GeoReport as psu

pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/LeicippusTesting/Analysis/'
edDataPath = 'C:/Dev/Github/ProteinDataFiles/ccp4_data/'
printPath = 'C:/Dev/Github/ProteinDataFiles/LeicippusTesting/Analysis/'

OrigPdb = '7a6a'
DenAdjPdb = '7a6a_DenAdj'
LapAdjPdb = '7a6a_LapAdj'
DenAdjPdb_samp = '7a6a_SampledDenAdj'
LapAdjPdb_samp = '7a6a_SampledLapAdj'

#This only needs to be done once, then we have the pdb files
#fac.getCsvFromCppResults_PdbFiles(printPath + '7a6a_SAMPLESFILE.csv',printPath + 'pdb7a6a_SampledDenAdj.ent', printPath + 'pdb7a6a_SampledLapAdj.ent')
#fac.getCsvFromCppResults_PdbFiles(printPath + '7a6a_ATOMSADJUSTEDFILE.csv',printPath + 'pdb7a6a_DenAdj.ent', printPath + 'pdb7a6a_LapAdj.ent')


georep = psu.GeoReport([OrigPdb,DenAdjPdb,DenAdjPdb_samp], pdbDataPath, edDataPath, printPath, ed=False, dssp=False)
#georepLap = psu.GeoReport([OrigPdb,LapAdjPdb,LapAdjPdb_samp], pdbDataPath, edDataPath, printPath, ed=False, dssp=False)
geos = ['C:O','C:N+1','N:N','O:O','C:C']
hueList = ['bfactor','rid','atomNo']
data = georep.getGeoemtryCsv(geos, hueList)
data['rid'] = data['rid'].astype(str)
data['ID'] = data['pdbCode'] + data['chain'] + data['rid'] + data['aa']
data.to_csv(printPath + "7a6a_data.csv", index=False)
dataOrig = data.query('pdbCode == "7a6a"')
dataAdj = data.query('pdbCode == "7a6a_denadj"')
dataSamp = data.query('pdbCode == "7a6a_sampleddenadj"')
dataOrig['ID'] = dataOrig['pdbCode'] + dataOrig['chain'] + dataOrig['rid'] + dataOrig['aa']
dataAdj['ID'] = dataAdj['pdbCode'] + dataAdj['chain'] + dataAdj['rid'] + dataAdj['aa']
dataAdj['N:N:Adj'] = dataAdj['N:N']
dataAdj['C:C:Adj'] = dataAdj['C:C']
dataAdj['O:O:Adj'] = dataAdj['O:O']
dataSamp['ID'] = dataSamp['pdbCode'] + dataSamp['chain'] + dataSamp['rid'] + dataSamp['aa']
dataSamp['N:N:Samp'] = dataSamp['N:N']
dataSamp['C:C:Samp'] = dataSamp['C:C']
dataSamp['O:O:Samp'] = dataSamp['O:O']
#make the pdb data in shape
dataOrigDiff = dataOrig[['ID','N:N','C:C','O:O']]
dataAdjDiff = dataAdj[['ID','N:N:Adj','C:C:Adj','O:O:Adj']]
dataSampDiff = dataSamp[['ID','N:N:Samp','C:C:Samp','O:O:Samp']]
#Join for diffs
dataDiff = dataOrigDiff.set_index('ID').join(dataAdjDiff.set_index('ID'))
dataDiff = dataDiff.set_index('ID').join(dataSampDiff.set_index('ID'))
dataDiff = dataDiff.dropna()
dataDiff.to_csv(printPath + "7a6a_diffs.csv", index=False)


#dataLap = georepLap.getGeoemtryCsv(geos, hueList)
georep.addHistogram(data=data,geoX='C:O', restrictions={'pdbCode':'7a6a'},hue='ID')
georep.addHistogram(data=data,geoX='C:O', restrictions={'pdbCode':'7a6a_denadj'},hue='ID')
georep.addHistogram(data=data,geoX='C:O', restrictions={'pdbCode':'7a6a_sampleddenadj'},hue='ID')

georep.addHistogram(data=data,geoX='C:N+1', restrictions={'pdbCode':'7a6a'},hue='ID')
georep.addHistogram(data=data,geoX='C:N+1', restrictions={'pdbCode':'7a6a_denadj'},hue='ID')
georep.addHistogram(data=data,geoX='C:N+1', restrictions={'pdbCode':'7a6a_sampleddenadj'},hue='ID')

georep.addHistogram(data=data,geoX='N:N', restrictions={'pdbCode':'7a6a'},hue='ID')
georep.addHistogram(data=data,geoX='N:N', restrictions={'pdbCode':'7a6a_denadj'},hue='ID')
georep.addHistogram(data=data,geoX='N:N', restrictions={'pdbCode':'7a6a_sampleddenadj'},hue='ID')

georep.addHistogram(data=data,geoX='C:C', restrictions={'pdbCode':'7a6a'},hue='ID')
georep.addHistogram(data=data,geoX='C:C', restrictions={'pdbCode':'7a6a_denadj'},hue='ID')
georep.addHistogram(data=data,geoX='C:C', restrictions={'pdbCode':'7a6a_sampleddenadj'},hue='ID')

georep.addHistogram(data=data,geoX='O:O', restrictions={'pdbCode':'7a6a'},hue='ID')
georep.addHistogram(data=data,geoX='O:O', restrictions={'pdbCode':'7a6a_denadj'},hue='ID')
georep.addHistogram(data=data,geoX='O:O', restrictions={'pdbCode':'7a6a_sampleddenadj'},hue='ID')

georep.addScatter(data=data,geoX='C:O',geoY='C:N+1',hue='pdbCode',categorical=True)
georep.addScatter(data=data,geoX='C:O',geoY='C:N+1',hue='pdbCode',categorical=True,exclusions={'pdbCode':'7a6a_sampleddenadj'})
georep.addScatter(data=data,geoX='C:O',geoY='C:N+1',hue='pdbCode',categorical=True,exclusions={'pdbCode':'7a6a_denadj'})

georep.printToHtml('Sampled',3,'Resampled_adj')

