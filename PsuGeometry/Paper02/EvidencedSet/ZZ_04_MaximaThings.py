import _Helpers as help

M_maxima_fake,M_maxima_real = False, False
M_top_20_report = False
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
    pdbs = help.getList('TOP20', 0, 'Fo_ADJ')
    pdbs.sort()
    reduce=False
    tag = '_realFo20'
    maxa.maximaCompareReal('Fo', pdbs, tag)
    #maxa.maximaCompareReal('QUINTIC', pdbs, tag)
    #maxa.maximaCompareReal('HEPTIC', pdbs, tag)
    #maxa.maximaCompareReal('NONIC', pdbs, tag)

if M_top_20_report:
    pdbs = help.getList('TOP20', 0, 'Fo_ADJ')
    realCsv, badRealCsv, occRealCsv = help.getMaximaDiffs('Fo', pdbs, False)
    print(occRealCsv)
    badAtoms = help.getBadList(realCsv, badRealCsv,occRealCsv, 0.5)
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
    data3 = help.getCsv('Fo' + '_ADJ', pdbs, geos, badAtoms, True, True, aa='ALL', includeCis=False, allAtoms=False,bFactorFactor=1.3, cutoff=0)
    #data5 = help.getCsv('QUINTIC' + '_ADJ', pdbs, geos,badAtoms, True, True, aa='ALL', includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=0)
    #data7 = help.getCsv('HEPTIC' + '_ADJ', pdbs, geos,badAtoms, True, True, aa='ALL', includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=0)
    #data9 = help.getCsv('NONIC' + '_ADJ', pdbs, geos,badAtoms, True, True, aa='ALL', includeCis=False, allAtoms=False, bFactorFactor=1.3, cutoff=0)

    import O_EvidenceSets as evidence
    evidence.evidenceReports('AscSets', ['Unrestricted','Restricted (Occ=1 Bff=1.3)','Restricted+= maxDiff=0.05','Adjusted Maxima 3 Fo'],dataPdbUn, dataPdb, dataPdbCut, data3, geoTrios, 'Evidential Geometry: Ascending Sets', perAA=False, tag='EH_SetsTopFo20_50')
    #evidence.evidenceReports('AscDegree', ['Adjusted Cubic','Adjusted Quintic','Adjusted Heptic','Adjusted Nonic'],data3, data5, data7, data9, geoTrios,'Evidential Geometry: Ascending Adjustments', perAA=False, tag='EH_DegTop20')

if F_compareCutAdjusted:
    import M_CompareMaxima as maxa
    pdbSet = 'Fo'
    geos = ['CA:C','C:O','C:N+1','N:CA']
    cutoff=0

    dataPdbCut = help.getCsv('RESTRICTED_CUT', [], [], [], False, False, aa='ALL', includeCis=False,allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
    dataAdj = help.getCsv(pdbSet + '_ADJ', [], [], [], False, False, aa='ALL', includeCis=False,allAtoms=False, bFactorFactor=1.3, cutoff=cutoff)
    combinedSet = help.comparisonDataSet(dataPdbCut,dataAdj,geos,'c:\dataFo.csv')
    maxa.compareAtomsPdbAdjusted(combinedSet, geos, pdbSet, 'Fo')
    setLower = combinedSet.query('RES <=0.7')
    setMiddle = combinedSet.query('RES <=0.85')
    setMiddle = setMiddle.query('RES > 0.7')
    setUpper = combinedSet.query('RES > 0.85')

    maxa.compareAtomsPdbAdjusted(setLower, geos, pdbSet, 'FoL')
    maxa.compareAtomsPdbAdjusted(setMiddle, geos, pdbSet, 'FoM')
    maxa.compareAtomsPdbAdjusted(setUpper, geos, pdbSet, 'FoU')
    #Take ut phenxix
    #combinedSet['SOFT'] = combinedSet['SOFTWARE'].str[:4]
    #combinedSet = combinedSet.query("SOFT != 'PHEN'")
