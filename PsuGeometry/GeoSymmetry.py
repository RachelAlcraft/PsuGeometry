
from PsuGeometry import GeoCalcs
from PsuGeometry import SymPoint



class GeoSymmetry:
    def __init__(self,peaks):
        self.distancePoints = {}
        #dataFrame(columns=('pdb_code', 'c', 'r', 's', 'x', 'y', 'z', '2FoFc', 'FoFc', 'Fo', 'Fc'))
        self.peaks = peaks

        # 1. Create distance matrix of all points and arrange into a distance group
        self.createDistances()


        # 2.l Put distances into groups of pairs

    def createDistances(self):
        cs = self.peaks['c']
        rs = self.peaks['r']
        ss = self.peaks['s']
        vs = self.peaks['Fc']
        count = 0
        for i in range(0,len(self.peaks)):
            for j in range(i+1, len(self.peaks)):
                ci,ri,si,vi = cs[i],rs[i],ss[i],vs[i]
                cj, rj, sj,vj = cs[j],rs[j],ss[j],vs[j]
                #cj,rj,sj = dfj['c'][0],dfj['r'][0],dfj['s'][0]
                dis = round(GeoCalcs.distance(ci,ri,si,cj,rj,sj),2)
                pi = SymPoint.SymPoint(ci,ri,si,vi)
                pj = SymPoint.SymPoint(cj, rj, sj, vj)
                if dis not in self.distancePoints:
                    self.distancePoints[dis] = []
                self.distancePoints[dis].append([pi,pj])
                if count%5000 == 0:
                    print('PSU symmetry',i,j,' ',count,'/',int((len(self.peaks)**2)/2))
                count = count+1

        for dis,pts in self.distancePoints.items():
            print(dis,len(pts),pts)



