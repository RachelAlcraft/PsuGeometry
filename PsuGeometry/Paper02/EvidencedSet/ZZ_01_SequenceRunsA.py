import gc

from matplotlib import pyplot as plt

import _Helpers as help


A_createGeos, B_merge = False,False
C_eh, D_mergeSets, E_compareSets = False, True,True


pdbSets = ['Fo']
cutoff = 0

if A_createGeos:
    import A_CreateGeosFileV2 as create2

    for pdbSet in pdbSets:
        pdbs = help.getList('70', cutoff, pdbSet + '_ADJ')
        #pdbs.remove('6euw')
        realCsv, badRealCsv, occRealCsv = help.getMaximaDiffs(pdbSet, pdbs, False)

        #create2.createGeosFile('RESTRICTED', pdbs, [], cutoff, '')
        #create2.createGeosFile('UNRESTRICTED', pdbs, [], cutoff, '')

        badAtoms = help.getBadList(realCsv, badRealCsv, occRealCsv,0.05)
        create2.createGeosFile(pdbSet + '_ADJ', pdbs,badAtoms, cutoff,'')
        create2.createGeosFile('RESTRICTED',pdbs, badAtoms, cutoff,'_CUTFo')


if B_merge:
    import B_MergeCsvs as merge
    #merge.mergeCsvs('RESTRICTED')
    merge.mergeCsvs('RESTRICTED_CUTFo')
    #merge.mergeCsvs('UNRESTRICTED')
    for pdbSet in pdbSets:
        merge.mergeCsvs(pdbSet + '_ADJ')

if C_eh:
    import C_EHCompare as compare
    compare.EHCompare('RESTRICTED')
    compare.EHCompare('RESTRICTED_CUTFo')
    compare.EHCompare('UNRESTRICTED')
    for pdbSet in pdbSets:
        compare.EHCompare(pdbSet + '_ADJ')

if D_mergeSets:
    pdbSets = ['Fo', 'CUBIC']
    pdbSetsM = ['UNRESTRICTED', 'RESTRICTED', 'RESTRICTED_CUT']
    for pdbSet in pdbSets:
        pdbSetsM.append(pdbSet + '_ADJ')


    import D_MergeSets as sets
    sets.mergeSets(pdbSetsM,'Fo_')

if E_compareSets:
    import E_ProduceSummaryReport as comp
    comp.compareSets('Fo_')



