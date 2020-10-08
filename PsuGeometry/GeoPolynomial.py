
import math

class GeoPolynomial:

    def __init__(self,vals,derivatives):
        self.vals = vals
        self.derivatives = derivatives
        self.diffCoeffs = self.buildDiffCoeffs(vals)
        self.polyCoeffs = self.buildPolyCoeffs(self.diffCoeffs)
        self.desc = self.stringDescription()

        ##TEST IT##
        #self.diffCoeffs = self.buildDiffCoeffs([1,4,7])
        #self.polyCoeffs = self.buildPolyCoeffs(self.diffCoeffs)
        #self.desc = self.stringDescription()
        #print(self.desc)
        #val1 = self.getValue(2,0)
        #val2 = self.getValue(2,1)
        #val3 = self.getValue(2,2)
        #print(val1,val2,val3)

        #print('poly',self.polyCoeffs)
        #print(self.getValue(1,0))
        #print(self.getValue(2,0))
        #print(self.getValue(0,0))

    def buildDiffCoeffs(self,vals):
        coeffs = []
        coeffs.append(vals[0])
        nextline = vals
        while len(nextline) > 1:
            newline = []
            for p in range(1,len(nextline)):
                lastpoint = nextline[p-1]
                thispoint = nextline[p]
                diffpoint = thispoint-lastpoint
                newline.append(diffpoint)
            coeffs.append(newline[0])
            nextline = newline
        for i in range(0,len(coeffs)):
            coeffs[i] = coeffs[i]/math.factorial(i)
        return coeffs

    def buildN_MinusCoeffs(self,degree):
        coeffs = [1]
        if degree > 0:
            for d in range(0,degree):
                rowval = (d+1)*-1
                newrow = coeffs[:] # careful to make a copy
                if len(coeffs)>1:
                    for c in range(0,len(coeffs)):
                        if c> 0:
                            cleft = coeffs[c-1]
                            cright = coeffs[c]
                            cnew = (cleft*rowval)+cright
                            newrow[c] = cnew
                cabove=coeffs[-1]
                cnew = cabove*rowval
                newrow.append(cnew)
                coeffs=newrow
        coeffs.reverse()
        return coeffs

    def buildPolyCoeffs(self,diffcoeffs):
        coeffs = []
        for i in range(0,len(diffcoeffs)):
            coeff = diffcoeffs[i]
            degree = i
            n_coeffs = self.buildN_MinusCoeffs(degree)
            for j in range(0,len(n_coeffs)):
                this_coeff = coeff * n_coeffs[j]
                if j >= len(coeffs):
                    coeffs.append(this_coeff)
                else:
                    coeffs[j] = coeffs[j] + this_coeff
        return coeffs

    def getValue(self,xval, differ):
        yval = 0
        start = differ # 0 means not differentiuated, 1 means first derivative etc
        for i in range(start,len(self.polyCoeffs)):
            origdegree = i
            degree = origdegree-differ
            coeff = self.polyCoeffs[i]*math.factorial(origdegree)/math.factorial(degree)
            newy = (xval**degree) * coeff
            yval = yval + newy
        return yval

    def stringDescription(self):
        desc = ''
        for i in range(0,len(self.polyCoeffs)):
            if desc != '':
                desc = desc + ' '
            desc = desc + str(self.polyCoeffs[i]) + 'x^' + str(i)
        return desc


