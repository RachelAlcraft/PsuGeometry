import gc

from matplotlib import pyplot as plt

import _Helpers as help

M_maxima_fake,M_maxima_real = False, False
M_top_20_report = False
A_createGeos, B_merge = False,False
C_eh, D_mergeSets, E_compareSets = False, False,False
F_compareCutAdjusted = True


if M_maxima_fake:
    import M_CompareMaxima as maxa
    pdbs = help.getList('TOP20', 0, 'CUBIC_FAKE_ADJ')
    #pdbs.remove('4rek') #forgot to reduce bfactor
    pdbs.sort()

    reduce=False
    tag = '_20'
    maxa.maximaCompareFake('CUBIC_FAKE', pdbs, tag,reduce)
    #maxa.maximaCompareFake('QUINTIC_FAKE', pdbs, tag,reduce)
    #maxa.maximaCompareFake('HEPTIC_FAKE', pdbs, tag,reduce)
    #maxa.maximaCompareFake('NONIC_FAKE', pdbs, tag,reduce)
    reduce = True
    tag = '_red_20'
    #maxa.maximaCompareFake('BOXES0', pdbs, tag, reduce)
    #maxa.maximaCompareFake('BOXES1', pdbs, tag, reduce)
    #maxa.maximaCompareFake('BOXES3', pdbs, tag, reduce)
    #maxa.maximaCompareFake('BOXES5', pdbs, tag, reduce)
    #maxa.maximaCompareFake('BOXES7', pdbs, tag, reduce)
    #maxa.maximaCompareFake('BOXES9', pdbs, tag, reduce)

if M_maxima_real:
    import M_CompareMaxima as maxa
    #pdbs = help.getList('70', 0, 'BOXES7_ADJ')
    pdbs = help.getList('TOP20', 0, 'CUBIC_ADJ')
    pdbs.sort()
    reduce=False
    tag = '_real20'
    maxa.maximaCompareReal('CUBIC', pdbs, tag)
    #maxa.maximaCompareReal('QUINTIC', pdbs, tag)
    #maxa.maximaCompareReal('HEPTIC', pdbs, tag)
    #maxa.maximaCompareReal('NONIC', pdbs, tag)

if M_top_20_report:
    pdbs = help.getList('TOP20', 0, 'CUBIC_ADJ')
    realCsv, badRealCsv, occRealCsv = help.getMaximaDiffs('CUBIC', pdbs, False)
    print(occRealCsv)
    badAtoms = help.getBadList(realCsv, badRealCsv,occRealCsv, 0.05)
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

    dataPdbUn = help.getCsv('PDB', pdbs, geos, [], True, True, aa='ALL', includeCis=False, allAtoms=True,bFactorFactor=-1, cutoff=0)
    dataPdb = help.getCsv('PDB', pdbs, geos, [], True, True, aa='ALL', includeCis=False, allAtoms=False,bFactorFactor=1.3, cutoff=0)
    dataPdbCut = help.getCsv('PDB', pdbs, geos, badAtoms, True, True, aa='ALL', includeCis=False, allAtoms=False,bFactorFactor=1.3, cutoff=0)
    data3 = help.getCsv('CUBIC' + '_ADJ', pdbs, geos, badAtoms, True, True, aa='ALL', includeCis=False, allAtoms=False,bFactorFactor=1.3, cutoff=0)
    data5 = help.getCsv('QUINTIC' + '_ADJ', pdbs, geos,badAtoms, True, True, aa='ALL', includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=0)
    data7 = help.getCsv('HEPTIC' + '_ADJ', pdbs, geos,badAtoms, True, True, aa='ALL', includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=0)
    data9 = help.getCsv('NONIC' + '_ADJ', pdbs, geos,badAtoms, True, True, aa='ALL', includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=0)

    import O_EvidenceSets as evidence
    evidence.evidenceReports('AscSets', ['Unrestricted','Restricted (Occ=1 Bff=1.3)','Restricted+= maxDiff=0.05','Adjusted Maxima Quintic'],dataPdbUn, dataPdb, dataPdbCut, data5, geoTrios, 'Evidential Geometry: Ascending Sets', perAA=False, tag='EH_SetsTop20')
    evidence.evidenceReports('AscDegree', ['Adjusted Cubic','Adjusted Quintic','Adjusted Heptic','Adjusted Nonic'],data3, data5, data7, data9, geoTrios,'Evidential Geometry: Ascending Adjustments', perAA=False, tag='EH_DegTop20')

pdbSets = ['CUBIC']
cutoff = 0

if A_createGeos:
    import A_CreateGeosFileV2 as create2

    for pdbSet in pdbSets:
        pdbs = help.getList('70', cutoff, pdbSet + '_ADJ')
        pdbs.remove('6euw')
        realCsv, badRealCsv, occRealCsv = help.getMaximaDiffs(pdbSet, pdbs, False)

        create2.createGeosFile('RESTRICTED', pdbs, [], cutoff, '')
        create2.createGeosFile('UNRESTRICTED', pdbs, [], cutoff, '')

        badAtoms = help.getBadList(realCsv, badRealCsv, occRealCsv,0.05)
        create2.createGeosFile(pdbSet + '_ADJ', pdbs,badAtoms, cutoff,'')
        create2.createGeosFile('RESTRICTED',pdbs, badAtoms, cutoff,'_CUT')


if B_merge:
    import B_MergeCsvs as merge
    merge.mergeCsvs('RESTRICTED')
    merge.mergeCsvs('RESTRICTED_CUT')
    merge.mergeCsvs('UNRESTRICTED')
    for pdbSet in pdbSets:
        merge.mergeCsvs(pdbSet + '_ADJ')

if C_eh:
    import C_EHCompare as compare
    compare.EHCompare('RESTRICTED')
    compare.EHCompare('RESTRICTED_CUT')
    compare.EHCompare('UNRESTRICTED')
    for pdbSet in pdbSets:
        compare.EHCompare(pdbSet + '_ADJ')

if D_mergeSets:
    pdbSetsM = ['UNRESTRICTED', 'RESTRICTED', 'RESTRICTED_CUT']
    for pdbSet in pdbSets:
        pdbSetsM.append(pdbSet + '_ADJ')


    import D_MergeSets as sets
    sets.mergeSets(pdbSetsM,'Boxes3_')

if E_compareSets:
    import E_ProduceSummaryReport as comp
    comp.compareSets('Boxes3_')

if F_compareCutAdjusted:
    import M_CompareMaxima as maxa
    pdbSet = 'CUBIC'
    geos = ['CA:C','C:O','C:N+1','N:CA']

    dataPdbCut = help.getCsv('RESTRICTED_CUT', [], [], [], False, False, aa='ALL', includeCis=False,allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
    dataAdj = help.getCsv(pdbSet + '_ADJ', [], [], [], False, False, aa='ALL', includeCis=False,allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
    combinedSet = help.comparisonDataSet(dataPdbCut,dataAdj,geos,'c:\data.csv')
    maxa.compareAtomsPdbAdjusted(combinedSet, geos, pdbSet, '')