'''
This class uses the pdb_eda library found here: https://pdb-eda.readthedocs.io/en/latest/index.html
Author: Rachel Alcraft
Date: 01/09/2020
Description:
Loads the matrices via the pdb_eda library and performs a simple normalisaiton (for the future make this configurable)
!!! Currently not the interpolated density as this is just a simple use of the library
'''
import pdb_eda
import numpy as np
import pandas as pd


class GeoDensity:

    def __init__(self,pdbCode,normalisation,pdbDataPath,edDataPath):
        # this defines the data we allow in the atom
        # Structure information
        self.pdbCode = pdbCode
        self.norm = normalisation
        pdb_eda.densityAnalysis.ccp4folder =edDataPath
        pdb_eda.densityAnalysis.pdbfolder = pdbDataPath
        self.analyser = pdb_eda.densityAnalysis.fromPDBid(pdbCode)
        self.factor = 1
        self.translation = 0
        try:
            alpha = self.analyser.densityObj.header.alpha
            self.valid = True
            if self.norm == 'fifty':
                # med = np.median(self.analyser.densityObj.density.ravel())
                minm = self.analyser.densityObj.density.min()
                maxm = self.analyser.densityObj.density.max()
                med = np.median(self.analyser.densityObj.density)
                print('PSU: density min=', minm, 'med=', med, 'max=', maxm)
                # med = np.mean(self.analyser.densityObj.density)

                self.translation = -1 * self.analyser.densityObj.density.min()
                self.factor =  50 / (med + self.translation)
                print('PSU: normalisation trans=',self.translation,'factor=',self.factor)
                print('PSU: normalisation min=', 0, 'med=', (med + self.translation) * self.factor, 'max=', (maxm + self.translation) * self.factor)

            print('PSU: created density for', self.pdbCode)
        except:
            print('PSU: !there is no density for', self.pdbCode)
            self.valid = False

    def getDensityXYZ(self,x,y,z): # this is not the intyerpolated density
        tFoFc = self.analyser.densityObj.getPointDensityFromXyz([x,y,z])
        tFoFc += self.translation
        tFoFc *= self.factor
        FoFc = self.analyser.diffDensityObj.getPointDensityFromXyz([x, y, z])
        FoFc *= self.factor
        Fo = tFoFc - FoFc
        Fc = tFoFc - 2*FoFc
        return [tFoFc,FoFc,Fo,Fc]

    def getDensityCRS(self,c,r,s): # this is not the intyerpolated density
        tFoFc = self.analyser.densityObj.getPointDensityFromCrs([c,r,s])
        tFoFc += self.translation
        tFoFc *= self.factor
        FoFc = self.analyser.diffDensityObj.getPointDensityFromCrs([c,r,s])
        FoFc *= self.factor
        Fo = tFoFc - FoFc
        Fc = tFoFc - 2*FoFc
        return [tFoFc,FoFc,Fo,Fc]


    def getPeaks(self,allPoints=False):
        if allPoints:
            print("PSU: Warning, the Density points function can take some minutes")
        else:
            print("PSU: Warning, the Density peaks function can take some minutes")
        matrix = self.analyser.densityObj.density
        maxMat = matrix.max()
        a, b, c = self.analyser.densityObj.density.shape
        print('\t\tPSU: Peaks=',a,'/',end=',')
        finalPeakList = []
        for i in range(0,a):
            peaked = True
            for j in range(0, b):
                peakList = []
                if allPoints:
                    peakList = self.getRowPoints(matrix, maxMat, i, j, -1)
                else:
                    peakList = self.getRowPeaks(matrix,i,j,-1)
                for peak in peakList:
                    usePoint = False
                    if allPoints:
                        usePoint = True
                    else:
                        usePoint = self.isPeak(matrix, peak[0], peak[1], peak[2], peak[3], )
                    if usePoint:
                        if peaked:
                            print(peak[0],end=',') # this is just to know it is working as this is very slow
                            peaked = False
                        c,r,s = peak[2],peak[1],peak[0] # the matrix coords seem to be inverted??
                        x,y,z = self.analyser.densityObj.header.crs2xyzCoord([c,r,s]) # convert to x,y,z coordinates to compare with the structure
                        #tfofc, fofc, fo, fc = self.getDensityXYZ(x,y,z)
                        tfofc, fofc, fo, fc = self.getDensityCRS(c,r,s)
                        finalPeakList.append([c,r,s,x,y,z,tfofc,fofc,fo,fc])

        densityData = pd.DataFrame(columns=('pdb_code', 'c', 'r', 's', 'x', 'y', 'z', 'peak2FoFc','peakFoFc','peakFo','peakFc'))
        print('', end='\n')
        print('\t\tPSU: Density complete, points=', len(finalPeakList), end='\n')
        for peak in finalPeakList:
            nextrow = len(densityData)
            densityData.loc[nextrow] = (
            self.pdbCode.upper(), peak[0], peak[1], peak[2], peak[3], peak[4], peak[5], peak[6],peak[7],peak[8],peak[9])  # switching ijk to crs


        return (densityData)

    def getRowPeaks(self,matrix,x,y,z):
        '''
        Gets the peaks for a row
        '''
        a,b,c = matrix.shape
        maxMat = matrix.max()

        divisor = 8
        if a < 120:
            divisor = 8
        elif a < 150:
            divisor = 8
        elif a < 190:
            divisor = 8

        #print(a,b,c)
        xRange = range(x,x+1)
        yRange = range(y,y+1)
        zRange = range(z,z+1)
        if x == -1:
            xRange = range(0,a)
        if y == -1:
            yRange = range(0,b)
        if z == -1:
            zRange = range(0,c)

        #print(xRange)
        #print(yRange)
        #print(zRange)

        peakList = []
        lastval = -1000
        lastCoordsVal = -1,-1,-1,0
        goingUp = False
        for i in xRange:
            for j in yRange:
                for k in zRange:
                    #print(i,j,k)
                    val = matrix[i,j,k]
                    if val > lastval:
                        goingUp = True
                    else:
                        if goingUp: # then we are now going back down
                            if abs(lastCoordsVal[3]) > maxMat/divisor:
                                peakList.append(lastCoordsVal)
                            goingUp = False

                    lastval = val
                    #lastCoordsVal = i,j,k,val
                    lastCoordsVal = i, j, k, val
        return (peakList)

    def getRowPoints(self,matrix,maxMat,x,y,z):
        '''
        Gets the peaks for a row
        '''
        a,b,c = matrix.shape
        maxMat = matrix.max()
        divisor = 4
        if a < 120:
            divisor = 8
        elif a < 150:
            divisor = 7
        elif a < 190:
            divisor = 5
        #print(a,b,c)
        xRange = range(x,x+1)
        yRange = range(y,y+1)
        zRange = range(z,z+1)
        if x == -1:
            xRange = range(0,a)
        if y == -1:
            yRange = range(0,b)
        if z == -1:
            zRange = range(0,c)
        peakList = []
        for i in xRange:
            for j in yRange:
                for k in zRange:
                    #print(i,j,k)
                    val = matrix[i,j,k]
                    if val >= maxMat/divisor:
                        lastCoordsVal = i, j, k, val
                        peakList.append(lastCoordsVal)
        return (peakList)

    def isPeak(self,matrix,x,y,z,val):
        a, b, c = matrix.shape
        # it is not a peak if all are lower, but same is ok (for now)
        xRange = range(x-1,x+2)
        yRange = range(y-1,y+2)
        zRange = range(z-1,z+2)
        isPeak = True
        for i in xRange:
            for j in yRange:
                for k in zRange:
                    if i>=0 and j>=0 and k>=0:
                        if i<a and j<b and k<c:
                            newval = matrix[i,j,k]
                            if newval > val:
                                isPeak = False


        return (isPeak)