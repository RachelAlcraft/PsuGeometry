
fileUnrestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataUnrestricted.csv'
fileRestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataRestricted.csv'
fileReduced = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataReduced.csv'
fileAdjusted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataAdjusted.csv'
fileMergedDiff = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csv'



def applyPdb(id):
    return id[:4]

import pandas as pd

dataCsvAdjusted = pd.read_csv(fileAdjusted)
dataCsvReduced = pd.read_csv(fileReduced)
dataCsvRestricted = pd.read_csv(fileRestricted)
dataCsvUnrestricted = pd.read_csv(fileUnrestricted)
dataCsvMerged = pd.read_csv(fileMergedDiff)

dataCsvAdjusted['pdbCode'] = dataCsvAdjusted.apply(lambda row: applyPdb(row['ID']), axis=1)
dataCsvReduced['pdbCode'] = dataCsvReduced.apply(lambda row: applyPdb(row['ID']), axis=1)
dataCsvRestricted['pdbCode'] = dataCsvRestricted.apply(lambda row: applyPdb(row['ID']), axis=1)
dataCsvUnrestricted['pdbCode'] = dataCsvUnrestricted.apply(lambda row: applyPdb(row['ID']), axis=1)
#dataCsvMerged['pdbCode'] = dataCsvMerged.apply(lambda row: applyPdb(row['id']), axis=1)

dataCsvAdjusted.to_csv(fileAdjusted, index=False)
dataCsvReduced.to_csv(fileReduced, index=False)
dataCsvRestricted.to_csv(fileRestricted, index=False)
dataCsvUnrestricted.to_csv(fileUnrestricted, index=False)
#dataCsvMerged.to_csv(fileMergedDiff, index=False)