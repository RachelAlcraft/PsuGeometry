
import pandas as pd
import math
import numpy as np

def getCsvFromCppResults_Slices(fileName):
    dataResults = {}
    #First load up as a big string
    text_file = open(fileName, "r")
    cpp_data = text_file.read()
    text_file.close()# close file

    startPos = cpp_data.find('BEGIN_SLICENUMS') + len('BEGIN_SLICENUMS')
    endPos = cpp_data.find('END_SLICENUMS')
    print(startPos,endPos)
    if endPos > startPos:
        numSlices = (cpp_data[startPos:endPos]).strip()
        print(numSlices)
        df = pd.DataFrame([numSlices])
        dataResults['SLICENUMS'] = df
        IDs = [['DENSITYSLICE','Density'],['RADIANTSLICE','Radiant'],['LAPLACIANSLICE','Laplacian']]

        for i in range(int(numSlices)):
            for ID,col in IDs:
                datalst = []
                ID_start = 'BEGIN_' + ID + '_' + str(i)
                startPos = cpp_data.find(ID_start) + len(ID_start)
                ID_end = 'END_' + ID + '_' + str(i)
                endPos = cpp_data.find(ID_end)
                data = (cpp_data[startPos:endPos]).strip()
                datalsta = data.split('\n')
                firstRow = True
                for row in datalsta:
                    if not firstRow:
                        lst = row.split(',')
                        datalst.append(lst)
                    firstRow=False
                df = pd.DataFrame(datalst, columns=['i', 'j', col])
                dataResults[ID + '_' + str(i)] = df


    return int(numSlices), dataResults

def getCsvFromCppResults_PdbFiles(fileName, denAdjusted, lapAdjusted):
    dataResults = {}
    #First load up as a big string
    text_file = open(fileName, "r")
    cpp_data = text_file.read()
    text_file.close()# close file

    startPos = cpp_data.find('BEGIN_DENSITYADJUSTED') + len('BEGIN_DENSITYADJUSTED')
    endPos = cpp_data.find('END_DENSITYADJUSTED')
    print(startPos,endPos)
    if endPos > startPos:
        denFile = (cpp_data[startPos:endPos]).strip()
        f = open(denAdjusted, "w")
        f.write(denFile)
        f.close()

    startPos = cpp_data.find('BEGIN_LAPLACIANADJUSTED') + len('BEGIN_LAPLACIANADJUSTED')
    endPos = cpp_data.find('END_LAPLACIANADJUSTED')
    print(startPos, endPos)
    if endPos > startPos:
        numSlices = (cpp_data[startPos:endPos]).strip()
        print(numSlices)
        lapFile = (cpp_data[startPos:endPos]).strip()
        f = open(lapAdjusted, "w")
        f.write(lapFile)
        f.close()


    return True


def DataFrameToMatrix(data, hue):
    real_len = len(data[hue].values)
    sq_len = int(math.sqrt(real_len))
    mtx = data[hue].values.reshape(int(sq_len), int(sq_len))
    npmtx = np.zeros((sq_len,sq_len))
    for i in range(sq_len):
        for j in range(sq_len):
            npmtx[i,j] = float(mtx[i,j])

    return npmtx
    #return np.array([[1, 2], [3, 4]])


