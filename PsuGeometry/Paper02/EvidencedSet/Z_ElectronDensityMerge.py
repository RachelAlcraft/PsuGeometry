

import H_MakeSlicesHtmlImage as slice

'''
*** ED Merging ***
'''


'''
ELECTRON DENSITY
'''
######################################################################################
#This can only be progressed to after the electron density has been calculated by Density Flight
######################################################################################
# H. Make html images fo the density slices

#slice.makeSlicesHtml('NCACO_B02_C_O', 'C:O Outliers, Set B02 (CA:C:O)')
#slice.makeSlicesHtml('NCACO_B02_CA_C', 'CA:C Outliers, Set B02 (CA:C:O)')

#slice.makeSlicesHtml('NCACO_B05_N_CA', 'N:CA Outliers, Set B05 (N:CA:C)')
#slice.makeSlicesHtml('NCACO_B05_TAU', 'TAU Outliers, Set B05 (N:CA:C)')
#slice.makeSlicesHtml('NCACO_B05_CA_C_O', 'CA:C:O Outliers, Set B05 (CA:C:O)')
#slice.makeSlicesHtml('RESTRICTED_TAU', 'TAU Outliers, Set Restricted (N:CA:C)')
#slice.makeSlicesHtml('UNRESTRICTED_TAU', 'TAU Outliers, Set Unrestricted (N:CA:C)')
#slice.makeSlicesHtml('NCACO_B04_TAU', 'TAU Outliers, Set B04 (N:CA:C)')
#slice.makeSlicesHtml('NCACO_B04_CACO', 'CA:C:O Outliers, Set B04 (CA:C:O)')
#slice.makeSlicesHtml('UNRESTRICTED_CACO', 'CA:C:O Outliers, Set Unrestricted (CA:C:O)')
#slice.makeSlicesHtml('NCACO_B02_CACO', 'CA:C:O Outliers, Set B02 (CA:C:O)','GoodOutliers')
#slice.makeSlicesHtml('NCACO_B03_NCA', 'N:CA Outliers, Set B03 (N:CA:C)','GoodOutliers')
#slice.makeSlicesHtml('B04_6q4g', 'Rejected from 6q4g, set B04 (CA:C:O)','Rejected')

'''
geoset.append([['CA', 'N', 'C'],['N:CA', 'CA:C', 'TAU']])
geoset.append([['C', 'CA', 'O'], ['C:O', 'CA:C:O']])
geoset.append([['CA', 'C-1', 'N'], ['C-1:N:CA']])
geoset.append([['C', 'O', 'N+1'], ['O:C:N+1']])
geoset.append([['C', 'CA', 'N+1'], ['C:N+1', 'CA:C:N+1']])
'''

#slice.makeSlicesHtml('NCACOSG_A0075_TAU', 'TAU Outliers, Set A0075 (N:CA:C)','GoodOutliers','NCACOSG_A0075_TAU')
slice.makeSlicesHtml('NCACO_ADJ01_ADJ_NCA', 'N:CA Outliers, Set ADJ01_ADJ (N:CA:C)','GoodOutliers','NCACO_ADJ01_ADJ_NCA')
#slice.makeSlicesHtml('NCACOSG_A0075_CAC', 'CA:C Outliers, Set A0075 (N:CA:C)','GoodOutliers','NCACOSG_A0075_CAC')

#slice.makeSlicesHtml('NCACOSG_A0075_CO', 'C:O Outliers, Set A0075 (CA:C:O)','GoodOutliers','NCACOSG_A0075_CO')
#slice.makeSlicesHtml('NCACOSG_A0075_CACO', 'CA:C:O Outliers, Set A0075 CA:C:O)','GoodOutliers','NCACOSG_A0075_CACO')

#slice.makeSlicesHtml('NCACOSG_A0075_C1NCA', 'C-1:N:CA (TAU-1) Outliers, Set A0075 (C-1:N:CA)','GoodOutliers','NCACOSG_A0075_C1NCA')

#slice.makeSlicesHtml('NCACOSG_A0075_OCN1', 'O:C:N+1 Outliers, Set A0075 (O:C:N+1)','GoodOutliers','NCACOSG_A0075_OCN1')

#slice.makeSlicesHtml('NCACOSG_A0075_CN1', 'C:N+1 Outliers, Set A0075 (CA:C:N+1)','GoodOutliers','NCACOSG_A0075_CN1')
#slice.makeSlicesHtml('NCACOSG_A0075_CACN1', 'CA:C:N+1 (TAU+1) Outliers, Set A0075 (CA:C:N+1)','GoodOutliers','NCACOSG_A0075_CACN1')




