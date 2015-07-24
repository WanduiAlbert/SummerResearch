import numpy as np
from gwpy.table.lsctables import SnglInspiralTable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#L1data = np.load('ER7L1triggers.npz', 'r')
#H1data = np.load('ER7H1triggers.npz', 'r')
import sys

fn = sys.argv[1]
gwdata = SnglInspiralTable.read(fn)
endtime =np.array(gwdata.getColumnByName('end_time')[:]) +\
    np.array(gwdata.getColumnByName('end_time_ns')[:]) * 1.0e-9
mass1 = np.array(gwdata.getColumnByName('mass1'))
mass2 = np.array(gwdata.getColumnByName('mass2'))
snr =  np.array(gwdata.getColumnByName('snr'))
chisq = np.array(gwdata.getColumnByName('chisq'))

mtotal = mass1 + mass2
eta = mass1 * mass2/ mtotal**2
mchirp = eta**(3./5)*mtotal

selection = snr >= 5.5 
plt.hist(snr[selection], bins=300)
ax = plt.gca()
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel('SNR')
#ax.set_xscale('log', nonposy='clip')
ax.set_ylabel('Number of events N')
ax.set_xlim(1, 5000)
plt.title('H1, BBH triggers, ER7')
plt.grid()
plt.savefig('H1_snr.png')


