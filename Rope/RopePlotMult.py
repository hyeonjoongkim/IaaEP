from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
pp   = PdfPages('RopePlotMult.pdf')
tmp1 = plt.figure(1)
plot = open('RopePlotMult-0.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
plt.plot( valx, valy, '-', label=r'CMS')
plot = open('RopePlotMult-1.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
plt.plot( valx, valy, '-', label=r'ALICE')
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,3))
plt.legend(frameon=False,loc='best')
plt.title(r'counts')
plt.xlabel(r'Multiplicity')
plt.ylabel(r'counts')
pp.savefig(tmp1,bbox_inches='tight')
plt.clf()
pp.close()
