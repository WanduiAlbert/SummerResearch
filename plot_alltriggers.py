import numpy as np
from gwpy.table.lsctables import SnglInspiralTable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as stats

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

# Make the boolean arrays for selecting the data that we need.
# Set a threshold on the SNR so that we can focus on the 
# fall off at lower SNRs.
min = 5.5
max = np.max(snr)
snr_sel = np.logical_and(snr >= min, snr <max)
mass_sel = np.vstack((np.logical_and(mchirp > 0.0, mchirp <= 5.0),\
    np.logical_and(mchirp > 5.0, mchirp <= 10.0), \
    np.logical_and(mchirp > 10.0, mchirp <= 15.0), \
    np.logical_and(mchirp > 15.0, mchirp <= 100.0)))
# Labels for the plots
labels = [r'$\mathcal{M}\ \leq\ 5M_{\odot}$', \
    r'$5M_{\odot}\ <\ \mathcal{M}\ \leq\ 10M_{\odot}$',\
   r'$10M_{\odot}\ <\ \mathcal{M}\ \leq\ 15M_{\odot}$',\
   r'$\mathcal{M}\ >\ 15M_{\odot}$']

data =[snr[np.logical_and(mass_sel[i,:], snr_sel)] for i in xrange(4)]
print len(data)
bins = np.logspace(np.log10(min), np.log10(max),num=50)
plt.figure(figsize=(12,10))
plt.hist(data, bins=bins, histtype='step', label=labels)
ax = plt.gca()
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel('SNR')
ax.set_xscale('log', nonposx='clip')
ax.set_ylabel('Number of events N')
ax.set_xlim(min, max)
plt.title('L1, BBH triggers, ER7')
plt.grid(True, which='both')
plt.legend()
plt.savefig('L1_snr_all_mchirp_bins_not_stacked.png')


