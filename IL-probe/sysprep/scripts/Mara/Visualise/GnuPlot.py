from os import system
from string import join

class GnuPlot:

    def __init__(self, path=""):
        self.path = path
        self.postfix = '.gnuplot'
        self.datafix = '.data'

    def plot(self, data, name=None):
        if name is None:
            name = 'test'
        dataname = self.path+name+self.datafix
        file = open(dataname,'w')
        c = len(data)
        l = len(data[0])
        for i in range(l):
            line = ""
            for j in range(c):
                line += str(data[j][i]) + " "
            file.write(line + '\n')
        file.close()
        plotname = self.path+name+self.postfix
        file = open(plotname,'w')
        plots = []
        for i in range(c-1):
            plots.append('"' + dataname + '" using ($1):($' + str(i+2) + ')')
        line = 'plot ' + join(plots,', ')
        file.write(line + '\npause -1 "Press <return>"\n')
        file.close()
        system('gnuplot ' + plotname + ' &')

