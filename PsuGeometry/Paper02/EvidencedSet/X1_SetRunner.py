

import A_CreateGeosFile as create
import B_MergeCsvs as merge
import C_EHCompare as compare
import I_ScatterReports as scatter
import F_CreateBadDensitySlicesCsv as bad
import G_CreateGoodOutlierSlicesCsv as good

'''
*** Paralleling ***
'''

pdbSets = ['NCACO_B02','NCACO_B025','NCACO_B03','NCACO_B04','NCACO_B05','RESTRICTED','UNRESTRICTED']
pdbSets = ['NCACO_B02','NCACO_B025','NCACO_B03']

for pdbSet in pdbSets:
    '''
    CREATION OF CSVS
    '''
    # A. Create the files of geoemtric measures - this takes a very long time
    #create.createGeosFile(pdbSet,0)
    # B. Merge the csv files together for comparative analaysis
    #merge.mergeCsvs(pdbSet)
    '''
    Engh&Huber Comparisons
    '''
    # C. Compare distributions with E&H values
    #compare.EHCompare(pdbSet)
    '''
    SCATTERS AND HISTOGRAMS
    '''
    # I. Scatters
    #scatter.scatterReports(pdbSet)


for pdbSet in pdbSets:
    '''
    ELECTRON DENSITY
    '''
    # F. Create density slices for the rejected density
    if pdbSet not in ['RESTRICTED', 'UNRESTRICTED']:
        bad.createBadDensitySlices(pdbSet, 'CA', 'N', 'C')
        bad.createBadDensitySlices(pdbSet, 'C', 'CA', 'O')
    # G. Create density slcies for the good outliers
    #good.createGoodDensitySlices(pdbSet, 'CA', 'N', 'C')
    #good.createGoodDensitySlices(pdbSet, 'C', 'CA', 'O')


