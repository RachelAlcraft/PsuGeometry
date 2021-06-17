

import _Helpers as help
import A_CreateGeosFile as create
import B_MergeCsvs as merge
import C_EHCompare as compare
import I_ScatterReports as scatter
import J_StatsCompare as stats
import K_StatsSummary as summary
import F_CreateBadDensitySlicesCsv as bad
import G_CreateGoodOutlierSlicesCsv as good

'''
*** Paralleling ***
'''

#pdbSets = ['NCACO_B02','NCACO_B025','NCACO_B03','NCACO_B04','NCACO_B05','RESTRICTED','UNRESTRICTED']
#pdbSets = ['NCACO_B02','NCACO_B025','NCACO_B03']
#pdbSets = ['NCACO_B02','NCACO_B025','NCACO_B03','NCACO_B04','NCACO_B05']
#pdbSets = ['NCACOSG_A01','NCACOSG_A0075','NCACOSG_A005','NCACOSG_A0025']
#pdbSets = ['NCACOSG_A0075']
pdbSets = ['CUBIC2']

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

    geos = ['C:N+1', 'N:N+1','TAU', 'PSI', 'PHI', 'C:O','O:C:N+1','CA:C:N+1','CA:C:O']
    pdbs = help.getList('TOP20',0)
    data = help.getCsv(pdbSet, pdbs,geos,True,True,aa='ALL',includeCis=False,allAtoms=False, bFactorFactor=1.3,cutoff=0)
    print(data)
    geoTrios = [['PHI', 'PSI', 'TAU'],
                 ['PSI', 'N:N+1', 'TAU'],
                 ['C:O', 'C:N+1', 'PHI'],
                ['C:O', 'C:N+1', 'PSI'],
                ['C:O', 'C:N+1', 'TAU'],
                ['PSI', 'C:N+1', 'TAU'],
                ['PSI', 'C:N+1', 'PHI'],
                ['PSI', 'C:O', 'TAU'],
                ['PSI', 'C:O', 'PHI'],
                ['PHI', 'C:N+1', 'TAU'],
                ['PHI', 'C:N+1', 'PSI'],
                ['PHI', 'C:O', 'TAU'],
                ['PHI', 'C:O', 'PSI'],
                ['C:N+1', 'C:O', 'O:C:N+1'],
                ['C:N+1', 'C:O', 'CA:C:O'],
                ['C:N+1', 'C:O', 'CA:C:N+1'],
                ['C:N+1', 'C:O', 'N:N+1'],
                ['C:N+1', 'N:N+1', 'TAU'],
                ]
    scatter.scatterReports(pdbSet,data,geoTrios,perAA =False,tag=pdbSet+'_tst')
    '''
    
    # J. Stats compare

    #stats.statsCompare(pdbSet,'RESTRICTED')
    #stats.statsCompare(pdbSet, 'UNRESTRICTED')

    # K. Stats summary
    #geos = ['CA:CB', 'CB:SG', 'N:CA:CB', 'CB:CA:C', 'CA:CB:SG', 'SG:{SG}']
    #summary.statsSummary(pdbSet, data, geos, '_DISULFIDE')
    #geos = ['N:CA', 'CA:C', 'C:O', 'C:N+1', 'TAU', 'CA:C:N+1', 'CA:C:O', 'O:C:N+1', 'C-1:N:CA']
    #data = help.getCsv(pdbSet, geos, False, True, 'ALL')
    #summary.statsSummary(pdbSet, data,geos,'EH')



for pdbSet in pdbSets:
    '''
    #ELECTRON DENSITY
    '''
    # F. Create density slices for the rejected density

    #if pdbSet not in ['RESTRICTED', 'UNRESTRICTED']:
    #    bad.createBadDensitySlices(pdbSet, 'CA', 'N', 'C')
    #    bad.createBadDensitySlices(pdbSet, 'C', 'CA', 'O')
    # G. Create density slcies for the good outliers
    geoset = []
    geoset.append([['CA', 'N', 'C'],['N:CA', 'CA:C', 'TAU']])
    geoset.append([['C', 'CA', 'O'], ['C:O', 'CA:C:O']])
    geoset.append([['CA', 'C-1', 'N'], ['C-1:N:CA']])
    geoset.append([['C', 'O', 'N+1'], ['O:C:N+1']])
    geoset.append([['C', 'CA', 'N+1'], ['C:N+1', 'CA:C:N+1']])
    good.createGoodDensitySlices(pdbSet, geoset)

'''




