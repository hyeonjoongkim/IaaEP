import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
pp   = PdfPages('RopePlot-allbins.pdf')
tmp1 = plt.figure(1)

data = [
	('RopePlot-0.dat','CMS, $N_{ch} < 35$',0.93),
	('RopePlot-1.dat','CMS, $35 < N_{ch} < 80$',0.2),
	('RopePlot-2.dat','CMS, $80 < N_{ch} < 105$',1),
	('RopePlot-3.dat','CMS $N_{ch} > 105$',1.1),
	('RopePlot-4.dat','ALICE, $N_{ch} < 35$',0.93),
	('RopePlot-5.dat','ALICE, $35 < N_{ch} < 80$',0.2),
	('RopePlot-6.dat','ALICE, $80 < N_{ch} < 105$',1),
	('RopePlot-7.dat','ALICE, $N_{ch} > 105$',1.1)
];
for i,(d1,d2) in enumerate(zip(data[0:4],data[4:8])):
	x = [None,None];
	y = [None,None];
	for j,d in enumerate((d1,d2)):
		with open(d[0],'r') as f:
			s = [line.split() for line in f];
			x[j] = np.array([float(t[0]) for t in s]);
			y[j] = np.array([float(t[1]) for t in s]);
	
	f1 = interpolate.interp1d(x[0],y[0]);
	f2 = interpolate.interp1d(x[1],y[1]);

	eps = 0.5;
	xs = np.linspace(1-eps,1+eps,1000);

	A = np.trapz(f1(xs),xs);
	B = np.trapz(f2(xs),xs);

	plt.plot(x[0],y[0]*d1[2],'-',label=d1[1]);
	plt.plot(x[1],y[1]*d2[2]*(A/B),':',label=d2[1],dashes=(1,2));

plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,3))
plt.legend(frameon=False,loc='best',prop={'size':6})
plt.title(r'CMS(2.0<|$\Delta\eta$|<4.0), ALICE(1.5<|$\Delta\eta$|<1.8)')
plt.xlabel(r'$\Delta\phi$ (rad)')
plt.ylabel(r'Abs. $(1/N_{\mathrm{N_Trigg}}) \mathrm{d}N / \mathrm{d}\Delta\phi$')
pp.savefig(tmp1,bbox_inches='tight')
plt.clf()
pp.close()
