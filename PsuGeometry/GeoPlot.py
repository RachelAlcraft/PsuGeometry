
import pandas as pd
import numpy as np
from matplotlib import cm
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import seaborn as sns
import io
import base64
import math

from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geopdb

class GeoPlot:
    def __init__(self,data,geoX,geoY='',title='',hue='bfactor',splitKey='',palette='viridis_r',
                 centre=False,vmin=0,vmax=0,operation='',newData=False,plot='scatter',categorical=False,
                 restrictions={},exclusions={},report=None,count=False):
        self.parent=report
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
        self.categorical = categorical
        self.numpy = []
        self.logged = False
        self.hasMatrix = False
        self.axX = 0,0
        self.axY = 0, 0
        self.differ=0
        self.restrictions=restrictions
        self.exclusions = exclusions
        if self.geoY == '' and plot not in 'surfaces':
            self.plot = 'histogram'
        self.count=count # only for histograms, probability or count
            #if self.hue=='bfactor':
            #    self.hue = 'pdbCode'
        if self.hue in 'aa,dssp,element,pdbCode':
            self.categorical=True

    #def getPlot(self,fig, ax):
    #    if self.plot == 'histogram':
    #        return self.plotHistogram(True,fig, ax)
    #    elif self.plot == 'scatter':
    #        return self.plotScatter(True,fig, ax)
    #    elif self.plot == 'probability':
    #        return self.plotProbability(True,fig, ax)

    def plotToAxes(self,fig, ax):
        if self.plot == 'histogram':
            return self.plotHistogram(fig, ax)
        elif self.plot == 'scatter' or self.plot=='contact':
            return self.plotScatter(fig, ax)
        elif self.plot == 'probability':
            return self.plotProbability(fig, ax)
        elif self.plot == 'surface':
            return self.plotSurface(fig, ax)
        elif self.plot == 'surfaces':
            return self.plotSurfaces(fig, ax)

    def plotSurface(self, fig, ax):
        afa = 1
        lw = 0.7
        if self.logged:
            afa = 1
        self.plotOneSurface(fig,ax,self.surface,afa,self.centre,self.palette,self.logged,lw)
        return ''

    def plotOneSurface(self, fig, ax,surface,afa,zero,palette,logged,lw):
        col='black'
        lvls=12
        if logged:
            col='black'
            x,y = surface.shape
            mind = 1000
            for i in range(0, x):
                for j in range(0, y):
                    mind = min(mind, surface[i,j])
            for i in range(0, x):
                for j in range(0, y):
                    val = (surface[i,j]-mind)+1
                    surface[i,j] = math.log(val)

        if zero:
            x, y = surface.shape
            mind = 1000
            maxd = -1000
            for i in range(0, x):
                for j in range(0, y):
                    mind = min(mind, surface[i, j])
                    maxd = max(maxd, surface[i, j])
            maxs = max(maxd,abs(mind))
            mins=-1*maxs
            image = plt.imshow(surface, cmap=palette, interpolation='nearest', origin='low', aspect='equal',vmin=mins,vmax=maxs,alpha=afa)
        else:
            image = plt.imshow(surface, cmap=palette, interpolation='nearest', origin='low', aspect='equal',alpha=afa)

        image = plt.contour(surface, colors=col, alpha=afa, linewidths=lw, levels=lvls)
        #ax.grid(False)
        #cbar = fig.colorbar(image, ax=ax)
        plt.axis('off')
        plt.title(self.title)
        return ''

    def plotSurfaces(self, fig, ax):
        afa = 0.95
        lw = 0.2
        for surface,palette,centre,logged in self.surface:
            self.plotOneSurface(fig,ax,surface,afa,centre,palette,logged,lw)
            afa = 0.55
            lw = 0.2
        return ''



    def plotHistogram(self,fig, ax):
        data = self.data.sort_values(by=self.geoX, ascending=True)
        title = self.title

        if self.operation == 'ABS':
            data = data[data[self.geoX] == abs(data[self.geoX])]
        elif self.operation == 'SQUARE':
            data = data[data[self.geoX] == data[self.geoX] ** 2]

        firstVal = data.head(1)[self.geoX].values[0]
        lastVal = data.tail(1)[self.geoX].values[0]
        firstHue = data.head(1)[self.hue].values[0]
        lastHue = data.tail(1)[self.hue].values[0]

        try:
            firstVal = round(firstVal, 2)
            lastVal = round(lastVal, 2)
        except:
            pass
        try:
            firstHue = round(firstHue, 2)
            lastHue = round(lastHue, 2)
        except:
            pass
        title += '\nFirst:' + self.hue + ' ' + str(firstHue) + '=' + str(firstVal)
        title += '\nLast:' + self.hue + ' ' + str(lastHue) + '=' + str(lastVal)

        # sns.distplot(data[xName], norm_hist=True, bins=50, kde=False)
        histCol = self.palette
        alpha=1
        bins = min(max(int(len(data[self.geoX])/6),10),50)
        density = not self.count
        if self.title == 'ghost':
            histCol = 'gainsboro'
            alpha=0.5
            plt.hist(data[self.geoX], EdgeColor='k', bins=bins,color=histCol,alpha=alpha,density=density,label='ghost')
            #sns.distplot(data[self.geoX], label='x', norm_hist=True, bins=50, kde=False,color='gainsboro')
        else:
            #if self.hue != '':
            #    splitList = data[self.hue].unique()
            #    for split in splitList:
            #        dfx = data[data[self.hue] == split]
            #        bins = max(int(len(dfx[self.geoX]) / 6),10)
            #        #plt.hist(dfx[self.geoX], EdgeColor='k', bins=bins, alpha=alpha, density=True)
            #        sns.distplot(dfx[self.geoX], label=split, norm_hist=True, bins=bins, kde=False,hist_kws=dict(alpha=0.5,EdgeColor='silver'))
            #    plt.legend()
            #else:
            #sns.distplot(data[self.geoX], label='', norm_hist=True, bins=bins, kde=False,hist_kws=dict(alpha=0.8,EdgeColor='silver'))
            plt.hist(data[self.geoX], EdgeColor='k', bins=bins, color=histCol, density=density,alpha=alpha)


        plt.title(title)
        plt.xlabel(self.geoX)

        if self.title != 'ghost':
            dfdesc = self.data[self.geoX].describe()
            rows = len(dfdesc.index)
            colsNames = list(dfdesc.index)
            html = "<table class='innertable'>\n"
            html += "<tr>\n"
            for r in range(0, rows):
                html += "<td>" + str(colsNames[r]) + "</td>\n"
            html += "</tr>\n"
            html += "<tr>"
            for r in range(0, rows):
                html += "<td>"
                try:
                    html += str(round(dfdesc[r], 2))
                except:
                    html += str(dfdesc[r])
                html += "</td>\n"

            html += "</tr>\n"
            html += "</table></p>\n"

            return html
        else:
            return ''

    def plotScatter(self,fig, ax):
        #fig, ax = plt.subplots()

        if self.categorical or self.hue == 'dssp':
            #blanksdata = self.data[self.data[self.hue] == '']
            #print(blanksdata)
            # it is possible for errors in dssp assignment in which case we call them X, but it will cover any errors not just dssp
            self.data.loc[self.data[self.hue] == '', self.hue] = 'X'

            gradients = {}
            gradsorig = self.data.sort_values(by=self.hue, ascending=True)[self.hue].unique()
            grads = self.getHueLists(self.hue,gradsorig)
            evenly_spaced_interval = np.linspace(0, 1, len(grads))

            try:
                sns.set_palette(sns.color_palette(self.palette, len(grads)))
                colors = [cm.get_cmap(self.palette)(x) for x in evenly_spaced_interval]
                i = 0
                for g in grads:
                    gradients[g] = colors[i]
                    i = i+1
                self.palette = gradients
            except:
                self.palette = self.palette

        if self.hue == 'resolution':
            self.data = self.data.sort_values(by=self.hue, ascending=False)
        elif self.plot == 'contact':
            self.data = self.data.sort_values(by='ridA', ascending=False)
        else:
            self.data = self.data.sort_values(by=self.hue, ascending=True)

        lw = 0.5
        alpha = 0.65
        ecol = 'grey'
        if self.palette == 'gist_gray_r':
            lw = 0  # this gives a crystollagraphic image look
            ecol = 'grey'
            alpha = 0.9
        if self.title=='ghost':
            alpha = 0.4
            self.palette='Greys'

        if self.centre:
            self.data[self.hue + '2'] = self.data[self.hue] ** 2
            data = self.data.sort_values(by=self.hue + '2', ascending=True)
            maxh = max(data[self.hue].max(), -1 * data[self.hue].min())
            minh = maxh * -1
            g = ax.scatter(data[self.geoX], data[self.geoY], c=data[self.hue], cmap=self.palette,
                           vmin=minh,vmax=maxh, edgecolor=ecol, alpha=alpha,linewidth=lw,s=20)
            cb = fig.colorbar(g)
            ax.set_xlabel(self.geoX)
            ax.set_ylabel(self.geoY)
            cb.set_label(self.hue)
        elif self.vmin < self.vmax:
            data = self.data.sort_values(by=self.hue, ascending=True)
            g = ax.scatter(data[self.geoX], data[self.geoY], c=data[self.hue], cmap=self.palette, vmin=self.vmin,
                           vmax=self.vmax, edgecolor=ecol, alpha=alpha,linewidth=lw,s=20)
            cb = fig.colorbar(g)
            ax.set_xlabel(self.geoX)
            ax.set_ylabel(self.geoY)
            cb.set_label(self.hue)

        elif self.plot == 'contact':
            alpha = 0.75
            self.data['distanceinv'] = 1/(self.data['distance'] ** 3)*4000
            if self.categorical == False:
                g = ax.scatter(self.data[self.geoX], self.data[self.geoY], c=self.data[self.hue],
                               cmap=self.palette,s=self.data['distanceinv'],edgecolor=ecol,alpha=alpha,linewidth=lw)
                cb = plt.colorbar(g)
                cb.set_label(self.hue)
            else:
                alpha = 0.65
                im = sns.scatterplot(x=self.geoX, y=self.geoY, hue=self.hue, data=self.data, alpha=alpha,legend='brief',
                                 palette=self.palette, size='distanceinv',edgecolor=ecol, linewidth=lw,vmax=3)
                #https://stackoverflow.com/questions/53437462/how-do-i-remove-an-attribute-from-the-legend-of-a-scatter-plot
                # EXTRACT CURRENT HANDLES AND LABELS
                h, l = ax.get_legend_handles_labels()
                # COLOR LEGEND (FIRST guess at size ITEMS) we don;t want to plot the distanceinc
                huelen = len(self.data.sort_values(by=self.hue, ascending=True)[self.hue].unique())+1
                col_lgd = plt.legend(h[:huelen], l[:huelen], loc='upper left',bbox_to_anchor=(1.05, 1), fancybox=True, shadow=True, ncol=1)
                plt.gca().add_artist(col_lgd)
                ax.set_xlabel('')
                ax.set_ylabel('')
        else:
            if self.categorical:
                im = sns.scatterplot(x=self.geoX, y=self.geoY, hue=self.hue, data=self.data, alpha=alpha,palette=self.palette
                                     ,edgecolor='aliceblue', linewidth=lw)
                if self.title!='ghost':
                    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  # Put the legend out of the figure
                else:
                    im.legend_.remove()
            else:
                g = ax.scatter(self.data[self.geoX], self.data[self.geoY], c=self.data[self.hue],
                               cmap=self.palette, edgecolor=ecol, alpha=alpha,linewidth=lw,s=20)
                cb = plt.colorbar(g)
                cb.set_label(self.hue)

        ax.set_xlabel(self.geoX)
        ax.set_ylabel(self.geoY)

        count = len(self.data.index)
        title = self.title
        if title == '':
            title += 'Count=' + str(count)
        else:
            title += '\nCount=' + str(count)

        plt.title(title)
        return ''

    def getPlotImage(self,fig, ax):
        #fig, ax = plt.subplots()
        img = io.BytesIO()
        fig.savefig(img, format='png', bbox_inches='tight')
        img.seek(0)
        encoded = base64.b64encode(img.getvalue())
        plt.close('all')
        return encoded

    def getAxes(self):
        xMin = min(self.data[self.geoX])
        xMax = max(self.data[self.geoX])
        yMin = min(self.data[self.geoY])
        yMax = max(self.data[self.geoY])
        return ([xMin,xMax,yMin,yMax])

    def plotProbability(self,fig, ax):

        # These shold be settings
        contours = 12
        kde = 0.10
        bins=50
        minX,maxX = self.axX[0],self.axX[1]
        minY,maxY = self.axY[0],self.axY[1]

        if minX == maxX:
            minX = min(self.data[self.geoX])
            maxX = max(self.data[self.geoX])
            minY = min(self.data[self.geoY])
            maxY = max(self.data[self.geoY])


        #fig, ax = plt.subplots()
        plt.axis([minX, maxX, minY, maxY])
        if self.hasMatrix:
            xgrid, ygrid, zgrid = self.numpy
        else:
            xgrid, ygrid, zgrid = self.kde2D_scipy(kde,[minX,maxX,minY,maxY], bins)

        xgrid = np.linspace(minX, maxX, bins)
        ygrid = np.linspace(minY,maxY, bins)

        ax.grid(True, which='major', axis='both', linestyle='-', color=(0.8, 0.8, 0.8), alpha=0.3)

        if self.centre:
            self.data[self.hue + '2'] = self.data[self.hue] ** 2
            self.data = self.data.sort_values(by=self.hue + '2', ascending=True)
            self.vmax = max(self.data[self.hue].max(), -1 * self.data[self.hue].min())
            self.vmin = self.vmax * -1

        alpha=1
        if self.title=='ghost':
            alpha = 0.4
            #self.palette = 'seismic'

        if self.vmin == self.vmax:
            im = plt.pcolormesh(xgrid, ygrid, zgrid, shading='gouraud', cmap=self.palette,alpha=alpha)
            cs = plt.contour(xgrid, ygrid, zgrid, contours, colors='0.7', linewidths=0.4,alpha=alpha)
        else:
            im = plt.pcolormesh(xgrid, ygrid, zgrid, shading='gouraud', cmap=self.palette, vmin=self.vmin, vmax=self.vmax,alpha=alpha)
            cs = plt.contour(xgrid, ygrid, zgrid, contours, colors='tab:purple', linewidths=0.05,alpha=alpha)

        cbar = fig.colorbar(im, ax=ax)
        cbar.remove()


        ax.set_xlabel(self.geoX)
        ax.set_ylabel(self.geoY)

        if self.hasMatrix:
            if self.title == '':
                title = 'Difference Image'
            else:
                title = self.title + '\nDifference Image'
        else:
            count = len(self.data.index)
            title = self.title
            if title == '':
                title += 'Count=' + str(count)
            else:
                title += '\nCount=' + str(count)

        plt.title(title)
        return ''

    def kde2D_scipy(self,bandwidth, axes, bins):

        xdata = self.data[self.geoX]
        ydata = self.data[self.geoY]
        data = np.vstack([xdata, ydata])
        xgrid = np.linspace(axes[0], axes[1], bins)
        ygrid = np.linspace(axes[2], axes[3], bins)
        Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
        grid_sized = np.vstack([Xgrid.ravel(), Ygrid.ravel()])
        # fit an array of size [Ndim, Nsamples]

        kde = gaussian_kde(data, bw_method=bandwidth)
        # evaluate on a regular grid
        Z = kde.evaluate(grid_sized)
        zgrid = Z.reshape(Xgrid.shape)
        return xgrid, ygrid, zgrid

    def getDifference(self,comparePlot):
        return ''

    def getOverlay(self,compareList):
        return ''

    def getHueLists(self,hue,huelist):
        if 'dssp' in hue: # could be be dsspA or dsspB too
            return ['-','B','E','G','H','I','S','T','X']
        else:
            return huelist

    def getNewData(self,hues=None):
        if self.plot == 'histogram':
            calcList = [self.geoX]
        else:
            calcList = [self.geoX, self.geoY]
        hueList = hues
        if hues == None:
            hueList = [self.hue]
        for rest in self.restrictions:
            if rest not in hueList:
                hueList.append(rest)

        dfs = []
        if self.parent != None:
            pdbmanager = geopdb.GeoPdbs(self.parent.pdbDataPath, self.parent.edDataPath, self.parent.ed, self.parent.dssp)
            for pdb in self.parent.pdbCodes:
                apdb = pdbmanager.getPdb(pdb)
                data = apdb.getGeoemtryCsv(calcList, hueList)
                dfs.append(data)
            self.data = pd.concat(dfs, ignore_index=True)
        else:
            print('PSU: cannot create data, pass in report parent to plots')

    def applyRestrictions(self):
        # now the data can be restricted as per the restrictions, which is a dictionary of restrictions, eg aa:'THR,PRO'
        if len(self.restrictions)>0 and not self.data.empty:
            dfs = []
            for rest in self.restrictions:
                allowed = self.restrictions[rest]
                allows = allowed.split(',')
                for all in allows:
                    data = self.data[self.data[rest] == all]
                    dfs.append(data)
                if self.title != '':
                    self.title += '\n'
                self.title += rest + '=' + allowed
            self.data = pd.concat(dfs, ignore_index=True)

    def applyExclusions(self):
        # now the data can be restricted as per the restrictions, which is a dictionary of restrictions, eg aa:'THR,PRO'
        if len(self.exclusions)>0 and not self.data.empty:
            dfs = []
            for exc in self.exclusions:
                notallowed = self.exclusions[exc]
                nallows = notallowed.split(',')
                for nall in nallows:
                    self.data = self.data[self.data[exc] != nall]
                if self.title != '':
                    self.title += '\n'
                self.title += exc + '!=' + notallowed

    def getMatrix(self):

        kde = 0.10
        bins = 50
        minX,maxX,minY,maxY = self.axX[0],self.axX[1],self.axY[0],self.axY[1]
        xgrid, ygrid, zgrid = self.kde2D_scipy(kde, [minX, maxX, minY, maxY], bins)
        return xgrid, ygrid, zgrid



class GeoOverlay:
    def __init__(self,plotA, plotB, title,report):
        self.title = title
        if title!='ghost':
            self.plotA = plotA
            self.plotB = plotB
        else:#In this case we have only the main plot, so we create the dummy plot
            self.plotB = plotA
            ghostReport = geor.GeoReport(['ghost'],report.pdbDataPath,report.edDataPath,report.outDataPath,report.ed,report.dssp)
            geoList = []
            geoList.append(self.plotB.geoX)
            if self.plotB.geoY != '':
                geoList.append(self.plotB.geoY)
            ghostdata = ghostReport.getGeoemtryCsv(geoList, ['pdbCode'])
            self.plotA = GeoPlot(ghostdata, self.plotB.geoX, geoY=self.plotB.geoY, title='ghost', hue='pdbCode', palette='Greys',plot=self.plotB.plot,operation=self.plotB.operation,report=report)

class GeoDifference:
    def __init__(self,dataA,dataB,geoX,geoY='',title='',restrictionsA={},restrictionsB={},exclusionsA={},exclusionsB={}, newData=False,palette='seismic',report=None):

        #huelist is all of the restrictions
        hues = []
        for hue in restrictionsA:
            if hue not in hues:
                hues.append(hue)


        self.plotA = GeoPlot(dataA,geoX,geoY=geoY,newData=newData,palette=palette,plot='probability',report=report)
        self.plotB = GeoPlot(dataB, geoX, geoY=geoY, newData=newData, palette=palette + '_r',plot='probability',report=report)

        self.plotA.restrictions = restrictionsA
        self.plotB.restrictions = restrictionsB

        self.plotA.exclusions = exclusionsA
        self.plotB.exclusions = exclusionsB

        if self.plotA.newData:
            self.plotA.getNewData(hues)

        if self.plotB.newData:
            self.plotB.getNewData(hues)

        titleA = self.plotA.title
        self.plotA.applyRestrictions() # these will be pplied twice, so don;t let the title get changed twice
        self.plotA.applyExclusions()
        self.plotA.title = titleA

        titleB = self.plotA.title
        self.plotB.applyRestrictions()
        self.plotB.applyExclusions()
        self.plotB.title = titleB

        self.plotA.newData = False
        self.plotB.newData = False

        axesA = self.plotA.getAxes()
        axesB = self.plotB.getAxes()
        axesAll = min(axesA[0],axesB[0]),max(axesA[1],axesB[1]),min(axesA[2],axesB[2]),max(axesA[3],axesB[3])
        self.plotA.axX = axesAll[0], axesAll[1]
        self.plotA.axY = axesAll[2], axesAll[3]
        self.plotB.axX = axesAll[0], axesAll[1]
        self.plotB.axY = axesAll[2], axesAll[3]

        arAA = self.plotA.getMatrix()
        arBB = self.plotB.getMatrix()

        arA = arAA[2]
        arB = arBB[2]
        arDiff = arA - arB
        minVal,maxVal=0,0
        for i in range(arA.shape[0]):
            for j in range(arA.shape[1]):
                maxVal = max(maxVal, arA[i, j])
                minVal = min(minVal, arA[i, j])

        for i in range(arB.shape[0]):
            for j in range(arB.shape[1]):
                maxVal = max(maxVal, arB[i, j])
                minVal = min(minVal, arB[i, j])

        maxVal = max(abs(maxVal), abs(minVal))
        minVal = 0 - maxVal

        self.plotA.vmin = minVal
        self.plotB.vmin = minVal
        self.plotA.vmax = maxVal
        self.plotB.vmax = maxVal


        self.plotDiff = GeoPlot(dataA, geoX, geoY=geoY, newData=False, palette=palette,vmin=minVal,vmax=maxVal,plot='probability',report=report)
        self.plotDiff.hasMatrix = True
        self.plotDiff.numpy = arB[0],arB[1],arDiff
        self.plotDiff.vmax = maxVal
        self.plotDiff.vmax = maxVal
        self.plotDiff.axX = axesAll[0], axesAll[1]
        self.plotDiff.axY = axesAll[2], axesAll[3]
