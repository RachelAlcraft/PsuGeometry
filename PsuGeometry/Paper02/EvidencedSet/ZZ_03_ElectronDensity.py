
runOutliers, runOutlierMerge = False,True
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
    #print(slices)
    good.createExtremeAdjustedSlices(pdbSet, slices,['C', 'CA', 'O'], 'C:O', dataFrame, '_O')
    good.createExtremeAdjustedSlices('RESTRICTED_CUT', slices, ['C', 'CA', 'O'], 'C:O', dataFrame, '_A')
    #the 2 files then need to be manually combined or copied into the overlay in density flight, with the tags changed

if runOutlierMerge:
    import H_MakeSlicesHtmlImage as slices
    slices.makeSlicesHtml('CO', 'AdjExtreme_CUBIC_ADJ_CO.csv','C:O Adjusted Extremes', 'COAdjExtremes','_3')