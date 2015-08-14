import numpy as np

import glob
import sys

from astropy.io import ascii

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

stats_dir = sys.argv[1]

channels = ["ALS-X_REFL_ERR_OUT_DQ",
"ASC-AS_B_RF36_I_YAW_OUT_DQ",
"ASC-AS_B_RF45_Q_PIT_OUT_DQ",
"ASC-AS_B_RF45_Q_YAW_OUT_DQ",
"ASC-Y_TR_B_NSUM_OUT_DQ",
"IMC-REFL_DC_OUT_DQ",
"PSL-ISS_AOM_DRIVER_MON_OUT_DQ"]

# get the files with the summary data written out
statsfiles = sorted(glob.glob(stats_dir + 'vetostats_[0-9][0-9].[0-9].txt'))

# this are the snr thresholds that were used in these calculations
snr_thresholds= np.arange(20,100,2)

# Now lets read out the data. We are interested in reading out the efficiency,
# deadtime and efficiency/deadtime columns

data = {}

for file in statsfiles:
  currdata = ascii.read(file, delimiter='|', header_start=0, data_start=1)
  for i in xrange(len(currdata)):
    channel = currdata['channel'][i]
    if channel not in data:
      data[channel] = [np.array(currdata['efficiency'][i]),\
          np.array(currdata['deadtime'][i]), np.array(currdata['efficiency/deadtime'][i])]
    else:
      data[channel] = map(lambda x, y: np.append(x, np.array(currdata[y][i])),\
          data[channel], ['efficiency', 'deadtime', 'efficiency/deadtime'])


print "All the data has been accurately read\n"

labels = ["%.1f" %_ for _ in snr_thresholds]

print "Starting to make the ROC curves"

for channel in data:
  plt.figure(figsize=(12,10))
  plt.plot(data[channel][1], data[channel][0], linestyle='-',\
      marker='o', color='b')
  for label, x,y in zip(labels[::4], data[channel][1][::4], data[channel][0][::4]):
    plt.annotate(label,
        xy = (x,y),
        xytext = (-20,20),
        textcoords = 'offset points', ha='right', va = 'bottom',
        bbox = dict(boxstyle = 'round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle='arc3,rad=0')
    )
  plt.title(r'ROC curve for channel %s' % channel)
  plt.xlabel('Deadtime [%]')
  plt.ylabel('Efficiency [%]')
  plt.grid(True, which='both')
  plt.savefig('%s_ROC.png' %channel)
  plt.close()

print "All done!\nMake efficiency/deadtime curves now"
for channel in data:
  plt.figure(figsize=(12,10))
  plt.plot(snr_thresholds, data[channel][2], linestyle='-',\
      marker='o', color='b')
  plt.title(r'%s' %channel)
  plt.grid(True, which='both')
  plt.xlabel('SNR Threshold')
  plt.ylabel('Efficiency/Deadtime')
  plt.savefig('%s_EffDt.png' %channel)
  plt.close()

print "Everything done!!!"
