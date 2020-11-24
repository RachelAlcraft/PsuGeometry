
from PsuGeometry import GeoTransformation as geot
import numpy as np

class GeoSpace:

    def getSquare(self,length,gap,central,linear,planar):
        print('PSU: get transformed square',central,linear,planar)
        # Create the transformation that would take us from whereever the atoms live to the origin.
        # We will then apply this transformation to all the coordinates in the cube
        transformation = geot.transformation(central, linear, planar)
        sq_x, sq_y, sq_z = self.generateEmptySquare(length,gap)
        count = 0
        for i in range(0, length):
            for j in range(0, length):
                    x = sq_x[i, j]
                    y = sq_y[i, j]
                    z = sq_z[i, j]
                    xt, yt, zt = transformation.applyTransformation([x, y, z])
                    sq_x[i, j] = xt
                    sq_y[i, j] = yt
                    sq_z[i, j] = zt

        return ([sq_x, sq_y, sq_z])

    def generateEmptySquare(self,length,gap):
        sq_x = np.zeros((length,length))
        sq_y = np.zeros((length,length))
        sq_z = np.zeros((length,length))
        sides = int((length-1)/2)

        for i in range(-1*sides,sides+1):
            for j in range(-1*sides,sides+1):
                    x = i * gap
                    y = j * gap
                    z = 0
                    sq_x[i+sides,j+sides] = x
                    sq_y[i+sides,j+sides] = y
                    sq_z[i+sides,j+sides] = z
        return ([sq_x,sq_y,sq_z])
