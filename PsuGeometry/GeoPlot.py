

class GeoQuery:
    def __init__(self,data,geoX,geoY='',title='',hue='2FoFc',splitKey='',palette='viridis_r',centre=False,vmin=0,vmax=0,operation='',newData=False,plot='scatter'):
        self.plot = plot
        self.data = data
        self.geoX = geoX
        self.geoY = geoY
        self.title=title
        self.hue = hue
        self.splitKey=splitKey
        self.palette=palette
        self.centre = centre
        self.vmin=vmin
        self.vmax=vmax
        self.operation=operation
        self.newData = newData
        if self.geoY == '':
            self.plot = 'histogram'

        def getPlot(self):
            if self.plot == 'histogram':
                return self.getHistogram()
            elif self.plot == 'scatter':
                return self.getScatter()
            elif self.plot == 'probability':
                return self.getProbability()




        def getScatter(self):
            return ''
        def getHistogram(self):
            return ''
        def getProbability(self):
            return ''
            












