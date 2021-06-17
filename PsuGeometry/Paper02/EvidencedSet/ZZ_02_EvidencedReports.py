import gc

from matplotlib import pyplot as plt

import _Helpers as help
import O_EvidenceSets as evidence
import C_EHCompareV2 as EH2
import J_StatsCompareV2 as stats2

runStats, runCompare, runScatters, runDifferenceHues, runAndReLoad = False,False,False,False,False

pdbSet = 'CUBIC'
cutoff = 0
reloadPDB = False
reloadCSV = False
extraTag = '3'

pdbs = help.getList('70', cutoff, pdbSet + '_ADJ')
#pdbs = ['1ucs','4zm7']
realCsv, badRealCsv, occRealCsv = help.getMaximaDiffs(pdbSet, pdbs, False)
badAtoms = help.getBadList(realCsv, badRealCsv, occRealCsv, 0.05)

if not reloadCSV:
    dataPdbUn = help.getCsv('UNRESTRICTED', pdbs, [], [], reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=True,bFactorFactor=-1, cutoff=cutoff)
    dataPdb = help.getCsv('RESTRICTED', pdbs, [], [], reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
    dataPdbCut = help.getCsv('RESTRICTED_CUT', pdbs, [],[],reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=False, bFactorFactor=1.3,cutoff=cutoff)
    dataAdj = help.getCsv(pdbSet + '_ADJ', pdbs, [], [], reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)

if runStats:
    #First set of reports are E&H histograms
    for pdbSet in [pdbSet]:
        geos = ['N:CA', 'CA:C', 'C:O', 'C:N+1', 'TAU', 'CA:C:N+1', 'CA:C:O', 'O:C:N+1', 'C-1:N:CA', 'CA-1:C-1:N:CA']
        geoTrios = [['N:CA'],
                    ['CA:C'],
                    ['C:O'],
                    ['C:N+1'],
                    ['TAU'],
                    ['CA:C:N+1'],
                    ['CA:C:O'],
                    ['O:C:N+1'],
                    ['C-1:N:CA'],
                    ]
        if reloadCSV:
            dataPdbUn = help.getCsv('PDB', pdbs, geos, [], reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=True, bFactorFactor=-1, cutoff=cutoff)
            dataPdb = help.getCsv('PDB', pdbs, geos, [], reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
            dataPdbCut = help.getCsv('PDB', pdbs, geos, badAtoms, reloadPDB, reloadCSV, aa='ALL', includeCis=False,allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
            dataAdj = help.getCsv(pdbSet + '_ADJ', pdbs, geos, badAtoms, reloadPDB, reloadCSV, aa='ALL',includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)

        evidence.evidenceReports(pdbSet,['Unrestricted','Restricted','+=Cut',pdbSet],dataPdbUn,dataPdb,dataPdbCut,dataAdj,geoTrios,'Evidential Geometry ' + pdbSet, perAA=False,tag='EH_70' + extraTag)

if runCompare:
    stats2.statsCompare('RESTRICTED_CUT', dataPdbCut,pdbSet + '_ADJ', dataAdj,True)

if runScatters:
    # Some scatters
    for pdbSet in [pdbSet]:
        geos = ['C:N+1', 'N:N+1', 'TAU', 'PSI', 'PHI', 'C:O', 'O:C:N+1', 'CA:C:N+1', 'CA:C:O']
        geoTrios = [['TAU'],
                    ['PHI', 'PSI', 'TAU'],
                    ['PSI', 'N:N+1', 'TAU'],
                    ['N:O-3', 'N:N+1', 'TAU'],
                    ['C:O-3', 'C:N+1', 'TAU'],
                    ['C:O', 'C:N+1', 'TAU'],
                    ['N:N+1', 'PHI', 'TAU'],
                    ['C-1:N:CA', 'CA:C:N+1', 'TAU'],
                    ]

        if reloadCSV:
            dataPdbUn = help.getCsv('PDB', pdbs, geos, [], reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=True, bFactorFactor=-1, cutoff=cutoff)
            dataPdb = help.getCsv('PDB', pdbs, geos, [], reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
            dataPdbCut = help.getCsv('PDB', pdbs, geos, badAtoms, reloadPDB, reloadCSV, aa='ALL', includeCis=False,allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
            dataAdj = help.getCsv(pdbSet + '_ADJ', pdbs, geos, badAtoms, reloadPDB, reloadCSV, aa='ALL',includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)


        evidence.evidenceReports(pdbSet,['Unrestricted','Restricted','+=Cut',pdbSet],dataPdbUn,dataPdb,dataPdbCut,dataAdj,geoTrios,'Evidential Geometry ' + pdbSet, perAA=False,tag='_scatters' + extraTag)

if runDifferenceHues:
    geos = ['C:N+1', 'N:N+1', 'TAU', 'PSI', 'PHI', 'C:O', 'O:C:N+1', 'CA:C:N+1', 'CA:C:O']
    geoTrios = [['TAU'],
                ['PHI', 'PSI', 'AvgMaximaDiffs'],
                ['PSI', 'N:N+1', 'AvgMaximaDiffs'],
                ['N:O-3', 'N:N+1', 'AvgMaximaDiffs'],
                ['C:O-3', 'C:N+1', 'AvgMaximaDiffs'],
                ['C:O', 'C:N+1', 'AvgMaximaDiffs'],
                ['N:N+1', 'PHI', 'AvgMaximaDiffs'],
                ['C-1:N:CA', 'CA:C:N+1', 'AvgMaximaDiffs'],
                ]

    if reloadCSV:
        dataPdbUn = help.getCsv('PDB', pdbs, geos, [], reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=True,bFactorFactor=-1, cutoff=cutoff)
        dataPdb = help.getCsv('PDB', pdbs, geos, [], reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=False,bFactorFactor=1.3, cutoff=cutoff)
        dataPdbCut = help.getCsv('PDB', pdbs, geos, badAtoms, reloadPDB, reloadCSV, aa='ALL', includeCis=False,allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
        dataAdj = help.getCsv(pdbSet + '_ADJ', pdbs, geos, badAtoms, reloadPDB, reloadCSV, aa='ALL', includeCis=False,allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)

    realCsv, badRealCsv, occRealCsv = help.getMaximaDiffs(pdbSet, pdbs, False)
    dataPdbUn = help.addMaximaDiffs(dataPdbUn,realCsv)
    dataPdb = help.addMaximaDiffs(dataPdb,realCsv)
    dataPdbCut = help.addMaximaDiffs(dataPdbCut,realCsv)
    dataAdj = help.addMaximaDiffs(dataAdj,realCsv)

    evidence.evidenceReports(pdbSet, ['Unrestricted','Restricted','+=Cut',pdbSet],dataPdbUn, dataPdb, dataPdbCut, dataAdj, geoTrios, 'Evidential Geometry ' + pdbSet, perAA=False, tag='_scattersDiffHue' + extraTag)

if runAndReLoad:
    reloadPDB = True
    reloadCSV = True
    geos = ['CB:O', 'N:O', 'N:C','TAU', 'PSI','PHI', 'N:N+1','CB:CA','N:CA:CB','CB:C:O']
    geoTrios = [
                ['PSI', 'N:N+1', 'TAU'],
                ['PSI', 'N:O', 'TAU'],
                ['CB:O', 'N:O', 'PSI'],
                ['CB:O', 'N:O', 'TAU'],
                ['N:N+1', 'N:O', 'PSI'],
                ['N:N+1', 'N:O', 'TAU'],
                ['N:N+1', 'N:O', 'PHI'],
                ['N:CA:CB', 'CB:C:O', 'PSI'],
                ['N:CA:CB', 'CB:C:O', 'TAU'],
                ['CB:CA', 'CB:O', 'PSI'],
                ['CB:CA', 'CB:O', 'TAU'],
                ['CB:C:O', 'CB:O', 'PSI'],
                ['CB:C:O', 'CB:O', 'TAU'],
                ['TAU', 'N:C', 'PSI'],
                ['TAU', 'N:C', 'PHI'],
                ['N:CA:CB', 'CB:O', 'PSI'],
                ['N:CA:CB', 'CB:O', 'TAU'],
                ['CB:C:O', 'CB:CA', 'PSI'],
                ['CB:C:O', 'CB:CA', 'TAU'],
                ['N:CA:CB', 'CB:CA', 'PSI'],
                ['N:CA:CB', 'CB:CA', 'TAU'],
                ]

    dataPdbUn = help.getCsv('PDB', pdbs, geos, [], reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=True,bFactorFactor=-1, cutoff=cutoff)
    dataPdb = help.getCsv('PDB', pdbs, geos, [], reloadPDB, reloadCSV, aa='ALL', includeCis=False, allAtoms=False,bFactorFactor=1.3, cutoff=cutoff)
    dataPdbCut = help.getCsv('PDB', pdbs, geos, badAtoms, reloadPDB, reloadCSV, aa='ALL', includeCis=False,allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
    dataAdj = help.getCsv(pdbSet + '_ADJ', pdbs, geos, badAtoms, reloadPDB, reloadCSV, aa='ALL', includeCis=False,allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
    #print('Unrestricted',dataPdbUn)
    #print('Restricted', dataPdb)
    #print('Cut Restricted', dataPdbCut)
    #print('Adjusted', dataAdj)

    evidence.evidenceReports(pdbSet,['Unrestricted','Restricted','+=Cut',pdbSet], dataPdbUn, dataPdb, dataPdbCut, dataAdj, geoTrios, 'Evidential Geometry ' + pdbSet, perAA=False, tag='_scattersSide' + extraTag)

