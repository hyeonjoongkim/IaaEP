from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
pp   = PdfPages('RopePlot.pdf')
tmp1 = plt.figure(1)
plot = open('RopePlot-0.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
plt.plot( valx, valy, '-', label=r'$N_ch < 35$')
plot = open('RopePlot-1.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
plt.plot( valx, valy, '-', label=r'$35< N_ch < 80$')
plot = open('RopePlot-2.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
plt.plot( valx, valy, '-', label=r'$80< N_ch < 105$')
plot = open('RopePlot-3.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
plt.plot( valx, valy, '-', label=r'$N_ch > 105$')
plot = open('RopePlot-4.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
plt.plot( valx, valy, '-', label=r'$N_ch < 35$')
plot = open('RopePlot-5.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
plt.plot( valx, valy, '-', label=r'$35< N_ch < 80$')
plot = open('RopePlot-6.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
plt.plot( valx, valy, '-', label=r'$80< N_ch < 105$')
plot = open('RopePlot-7.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]
plt.plot( valx, valy, '-', label=r'$N_ch > 105$')
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,3))
plt.legend(frameon=False,loc='best')
plt.title(r'dN/d$\Delta\phi$')
plt.xlabel(r'$\Delta\phi$ (rad)')
plt.ylabel(r'$(1/N_{\mathrm{N_Trigg}}) \mathrm{d}N / \mathrm{d}\Delta\phi$')
pp.savefig(tmp1,bbox_inches='tight')
plt.clf()
pp.close()
