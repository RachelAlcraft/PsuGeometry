
import D_MergeSets as sets
import E_ProduceSummaryReport as comp



'''
*** Merging ***
'''

pdbSets = ['NCACO_B02','NCACO_B025','NCACO_B03','NCACO_B04','NCACO_B05','RESTRICTED','UNRESTRICTED']

# D. Merge sets for comparison across sets
sets.mergeSets(pdbSets)

# E. Produce a summary report
comp.compareSets([])


