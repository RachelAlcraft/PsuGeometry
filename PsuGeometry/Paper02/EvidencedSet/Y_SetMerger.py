





'''
*** Merging ***
'''

#pdbSets = ['NCACO_B02','NCACO_B025','NCACO_B03','NCACO_B04','NCACO_B05','RESTRICTED','UNRESTRICTED']
pdbSets = ['NCACOSG_A0025','NCACOSG_A005','NCACOSG_A0075', 'NCACOSG_A01','RESTRICTED','UNRESTRICTED']
pdbSets = ['NCACOSG_A0075', 'NCACO_ADJ01','NCACO_ADJ01_ADJ','RESTRICTED','UNRESTRICTED']

# D. Merge sets for comparison across sets
import D_MergeSets as sets
sets.mergeSets(pdbSets,'AngAdj_')

# E. Produce a summary report
import E_ProduceSummaryReport as comp
comp.compareSets('AngAdj_')


