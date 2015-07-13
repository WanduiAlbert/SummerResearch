import numpy as np
from gwpy.table.lsctables import SnglInspiralTable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

L1data = np.load('ER7L1triggers.npz', 'r')
H1data = np.load('ER7H1triggers.npz', 'r')

endtime =np.array(L1data['arr_0'])
mass1 = np.array(L1data['arr_1'])
mass2 = np.array(L1data['arr_2'])
snr =  np.array(L1data['arr_3'])
chisq = np.array(L1data['arr_4'])

mtotal = mass1 + mass2
eta = mass1 * mass2/ mtotal**2
mchirp = eta**(3./5)*mtotal

selection = snr >= 5.5
plt.hist(snr[selection], bins=50)
plt.xlabel('SNR')
plt.ylabel('Number of events N')
plt.title('L1, BBH triggers, ER7')
plt.grid()
plt.savefig('L1_snr.png')

plt.figure()
plt.hist(mchirp[selection], bins=50)
plt.xlabel('Chirp Mass $\mathcal{M}$ $(M_{\odot})$')
plt.ylabel('Number of events N')
plt.title('L1, BBH triggers, ER7')
plt.grid()
plt.savefig('L1_mchirp.png')

plt.figure()
plt.scatter(mchirp[selection], snr[selection])
plt.xlabel('Chirp Mass $\mathcal{M}$ $(M_{\odot})$')
plt.ylabel('SNR')
plt.grid()
plt.title('L1, BBH triggers, ER7')
plt.savefig('L1_snr_vs_mchirp.png')
