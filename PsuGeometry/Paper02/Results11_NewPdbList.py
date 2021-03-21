# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import random
import pandas as pd
'''
TAU correlations
'''
###############################################################################################
myWindowsLaptop = True #only works on windows
bfactorFactor = 1.3
keepDisordered = False

geoList = ['N:N+1','TAU','PSI','PHI','N:C','CA:C','C:O','N:CA']
geoListTau = ['TAU']
hueList = ['aa', 'rid', 'bfactor','pdbCode']
aas = ['GLY']

runs = []
#runs.append('Better')
runs.append('Good')
#runs.append('Bad')
#runs.append('Orig')


pdbReDataPath = 'F:/Code/ProteinDataFiles/pdb_out/good/'
pdbdata = pd.read_csv('../PdbLists/Pdbs_Under1.csv') # This is a list of pdbs <= 1.1A non homologous to 90%

pdbListIn = pdbdata['PDB'].tolist()[0:]
pdbList = []
for pdb in pdbListIn:
    import os.path
    if os.path.isfile((pdbReDataPath + 'pdb' + pdb + '.ent').lower()):
        pdbList.append(pdb.lower())
    else:
        print('No file:',(pdbReDataPath + 'pdb' + pdb + '.ent').lower())

print(pdbList)

for dataSet in runs:
    ###################################################################################
    tag=dataSet
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'
    if dataSet == 'Good':
        pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_out/good/'
    elif dataSet == 'Better':
        pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_out/better/'
    elif dataSet == 'Bad':
        pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_out/bad/'



    #pdbList = ['6e6o']
    print('Creating ordered report')
    pdbmanager = geopdb.GeoPdbs(pdbDataPath, edDataPath, False,False,False)
    georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)


    print('Getting dataframes of geometry -- 1')
    dataA = georep.getGeoemtryCsv(geoList, hueList,bfactorFactor)
    dataATau = georep.getGeoemtryCsv(geoListTau, hueList, bfactorFactor)
    if dataATau.shape[0] > 0:
        #print(dataATau)

        for aa in aas:
            print(aa)
            sql = 'aa == "' + aa + '"'
            if dataA.shape[0] > 0:
                dataAaa = dataA.query(sql)
            else:
                dataAaa = dataA

            dataAaaTau = dataATau.query(sql)
            #print('TAU',dataAaaTau)

            if dataAaaTau.shape[0] > 0:

                if dataAaa.shape[0] > 0:
                    dataAmin = dataAaa.query('TAU > 100')
                    dataAmin = dataAmin.query('TAU < 125')


                dataAminTau = dataAaaTau.query('TAU > 100')
                dataAminTau = dataAminTau.query('TAU < 125')

                if dataAaa.shape[0] > 0:
                    georep.addScatter(data=dataAaa, geoX='PHI', geoY='PSI', hue='TAU', title='Ordered factor=1.3 ' + aa,palette='jet', sort='NON')
                georep.addHistogram(data=dataAaaTau, geoX='TAU', title='Ordered factor=1.3 ' + aa)
                if dataAaa.shape[0] > 0:
                    georep.addScatter(data=dataAaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='Ordered factor=1.3 ' + aa, palette='jet', sort='NON')

                if dataAaa.shape[0] > 0:
                    georep.addScatter(data=dataAmin, geoX='PHI', geoY='PSI', hue='TAU', title='Ordered factor=1.3 ' + aa,palette='jet', sort='NON')
                georep.addHistogram(data=dataAminTau, geoX='TAU', title='Ordered factor=1.3 ' + aa)
                if dataAaa.shape[0] > 0:
                    georep.addScatter(data=dataAmin, geoX='PSI', geoY='N:N+1', hue='TAU', title='Ordered factor=1.3 ' + aa, palette='jet', sort='NON')
                    georep.addScatter(data=dataAmin, geoX='N:C', geoY='CA:C', hue='C:O', title='' + aa, palette='jet', sort='NON')
                    georep.addScatter(data=dataAmin, geoX='N:CA', geoY='CA:C', hue='C:O', title='' + aa, palette='jet', sort='NON')

                    georep.addScatter(data=dataAmin, geoX='CA:C', geoY='C:O', hue='N:C', title='' + aa, palette='jet', sort='NON')
                    georep.addScatter(data=dataAmin, geoX='CA:C', geoY='C:O', hue='TAU', title='' + aa, palette='jet', sort='NON')

                    georep.addScatter(data=dataAmin, geoX='C:O', geoY='N:C', hue='CA:C', title='' + aa, palette='jet', sort='NON')
                    georep.addScatter(data=dataAmin, geoX='C:O', geoY='N:CA', hue='CA:C', title='' + aa, palette='jet',sort='NON')

                    georep.addScatter(data=dataAmin, geoX='N:C', geoY='CA:C', hue='TAU', title='' + aa, palette='jet', sort='NON')
                    georep.addScatter(data=dataAmin, geoX='N:CA', geoY='CA:C', hue='TAU', title='' + aa, palette='jet',sort='NON')


                    georep.addScatter(data=dataAmin, geoX='C:O', geoY='N:C', hue='TAU', title='' + aa, palette='jet', sort='NON')
                    georep.addScatter(data=dataAmin, geoX='C:O', geoY='N:CA', hue='TAU', title='' + aa, palette='jet', sort='NON')


                print('Creating reports')

                dicPdbs = []
                dicPdbsStd = []
                #now look at each pdb
                for pdb in pdbList:
                    sql = 'pdbCode == "' + pdb.lower() + '"'
                    dataPdb = dataAaaTau.query(sql)
                    dfdesc = dataPdb['TAU'].describe()
                    #print(dfdesc)
                    if dfdesc['count'] > 0:
                        dicPdb = {}
                        dicPdb['pdbCode'] = pdb
                        dicPdb['TAU_MEAN'] = dfdesc['mean']
                        dicPdb['TAU_MEDIAN'] = dfdesc['50%']
                        dicPdb['TAU_STD'] = dfdesc['std']
                        dicPdb['COUNT'] = dfdesc['count']
                        dicPdbs.append(dicPdb)

                    if dfdesc['count'] > 1:
                        dicStd = {}
                        dicStd['TAU_STD'] = 0
                        dicStd['pdbCode'] = pdb
                        dicPdbsStd.append(dicPdb)

                dataFramePdbs = pd.DataFrame.from_dict(dicPdbs)
                dataFramePdbsStd = pd.DataFrame.from_dict(dicPdbsStd)
                #print(dataFramePdbs)
                #print(dataFramePdbsStd)

                georep.addHistogram(data=dataFramePdbs, geoX='COUNT', title='Count per pdb ' + aa)
                georep.addHistogram(data=dataFramePdbs, geoX='TAU_MEDIAN', title='Median per pdb ' + aa)

                georep.addHistogram(data=dataFramePdbs, geoX='TAU_MEAN', title='Mean per pdb ' + aa)
                georep.addHistogram(data=dataFramePdbsStd, geoX='TAU_STD', title='Std per pdb ' + aa)

                georep.printToHtml('Results 11. Tau Plots\nPdbs=' + str(len(pdbList)) + '\nWith bFactorFactor of ' + str(bfactorFactor), 3, 'Results11_' + tag + aa)
                dataAaa.to_csv(printPath + "Results11" + tag + ".csv", index=False)
                dataAaaTau.to_csv(printPath + "Results11" + tag + "Tau.csv", index=False)
                dataFramePdbs.to_csv(printPath + "Results11" + tag + "_statsPdb.csv", index=False)

        else:
            print("No data for", aa)

    else:
        print("No data for",pdbDataPath)

    dataFramePdbs = None
    dataFramePdbsStd = None
    dataAaa = None
    dataAmin = None
    dataAaaTau = None
    dataAminTau = None
    georep = None
    pdbmanager.clear()





