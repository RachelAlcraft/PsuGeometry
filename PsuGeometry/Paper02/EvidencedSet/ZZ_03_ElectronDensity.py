
runOutliers, runOutlierMerge, runManualReport = False,False,True
if runOutliers:
    '''
    #ELECTRON DENSITY
    '''
    # G. Create density slices for the good outliers or extremes
    pdbSet = 'CUBIC_ADJ'
    import G_CreateGoodOutlierSlicesCsv as good
    geoset = []
    geoset.append([['C', 'CA', 'O'], ['C:O', 'CA:C:O']])
    slices, dataFrame = good.getExtremeAdjustedResidues(pdbSet, 'C:O',[1.173,1.317])
    good.produceExtremeValFile('RESTRICTED_CUT',slices,'C:O')
    #print(slices)
    good.createExtremeAdjustedSlices(pdbSet, slices,['C', 'CA', 'O'], 'C:O', dataFrame, '_O')
    good.createExtremeAdjustedSlices('RESTRICTED_CUT', slices, ['C', 'CA', 'O'], 'C:O', dataFrame, '_A')
    #the 2 files then need to be manually combined or copied into the overlay in density flight, with the tags changed

if runOutlierMerge:
    import H_MakeSlicesHtmlImage as slices
    slices.makeSlicesHtml('CO', 'OrigExtreme_RESTRICTED_CUT_CO.csv','AdjExtreme_CUBIC_ADJ_CO.csv','C:O Adjusted Extremes', 'COAdjExtremes','3')

if runManualReport:
    import H_MakeSlicesHtmlImage as slices
    dir = "C:/Dev/Github/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/SlicesI/"
    pdb, rid = '3w5h', 'x'
    pdb, rid = '6s2m', 'A11'
    lineRuns =[]
    lineRuns.append([pdb + ' ' + rid + ' Positions\n2Fo-Fc','cubehelix_r',pdb+'/Orig_2FoFc.csv',pdb+'/Orig_pos.csv'])
    lineRuns.append([pdb + ' ' + rid + ' Adjusted Positions\n2Fo-Fc', 'cubehelix_r', pdb+'/Adj_2FoFc.csv', pdb+'/Adj_pos.csv'])

    lineRuns.append(['Simulated Originial IAM', 'cubehelix_r', pdb+'/Orig_IAM.csv', pdb+'/Orig_pos.csv'])
    lineRuns.append(['Simulated Adjusted IAM', 'cubehelix_r', pdb+'/Adj_IAM.csv', pdb+'/Adj_pos.csv'])

    lineRuns.append(['Simulated Originial BEM', 'cubehelix_r', pdb+'/Orig_BE.csv', pdb+'/Orig_pos.csv'])
    lineRuns.append(['Simulated Adjusted BEM', 'cubehelix_r', pdb+'/Adj_BE.csv', pdb+'/Adj_pos.csv'])

    lineRuns.append(['Positions pdb vs adjusted', 'jet_r', pdb + '/Orig_pos.csv', pdb + '/Adj_pos.csv'])

    #Distance=1.227Å
    #lineRuns.append(['3w5h Fo Adjusted', 'cubehelix_r', 'Adj_Fo.csv', 'Adj_pos.csv'])
    #lineRuns.append(['3w5h Fc Adjusted', 'cubehelix_r', 'Adj_Fc.csv', 'Adj_pos.csv'])
    # lineRuns.append(['3w5h Fo', 'cubehelix_r', '3w5h/Orig_Fo.csv', '3w5h/Orig_pos.csv'])
    # lineRuns.append(['3w5h Fc', 'cubehelix_r', '3w5h/Orig_Fc.csv', '3w5h/Orig_pos.csv'])
    #lineRuns.append(['Simulated Adjusted BE 4H 0.2', 'cubehelix_r', 'Adj_BE2.csv', 'Adj_pos.csv'])
    # lineRuns.append(['Simulated Originial BE 4H 0.2', 'cubehelix_r', 'Orig_BE2.csv', 'Orig_pos.csv'])

    slices.makeSlicesHtmlFromValues(pdb + ' adjusted and original<br> EH 1991=1.231Å<br>Adj=1.241Å',dir,lineRuns,2,pdb)