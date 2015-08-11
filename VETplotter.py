# Import statements

import sys

import numpy as np

from gwpy.segments import Segment, SegmentList, DataQualityFlag
from gwpy.table.lsctables import SnglBurstTable, SnglInspiralTable

from glue.lal import Cache

from  gwpy.time import LIGOTimeGPS

from gwpy.toolkits.vet import (get_triggers,get_segments,get_metric)

import glob, os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Save all the channels in a python script
channels = ["ASC-AS_B_RF36_I_YAW_OUT_DQ",
"PEM-EY_MAG_EBAY_SUSRACK_QUAD_SUM_DQ",
"LSC-POP_A_LF_OUT_DQ",
"ASC-AS_B_RF45_Q_YAW_OUT_DQ",
"OMC-ASC_POS_Y_OUT_DQ",
"PEM-EX_TILT_VEA_FLOOR_T_DQ",
"PEM-EY_TILT_VEA_FLOOR_T_DQ",
"OMC-PZT2_MON_AC_OUT_DQ",
"LSC-SRCL_OUT_DQ",
"IMC-REFL_DC_OUT_DQ",
"ALS-X_REFL_ERR_OUT_DQ",
"PEM-CS_MAG_EBAY_SUSRACK_Y_DQ",
"ASC-Y_TR_B_NSUM_OUT_DQ",
"ASC-AS_B_RF45_Q_PIT_OUT_DQ",
"SUS-OMC_M1_ISIWIT_T_DQ",
"PSL-ISS_AOM_DRIVER_MON_OUT_DQ",
"LSC-PRCL_OUT_DQ"]

ifo = sys.argv[1]
bbhdir = sys.argv[2]
bbhfile= glob.glob(os.path.join(bbhdir, ifo+'*.xml.gz'))[0]
omiccachedir = sys.argv[3]

# Read in the segment file
segments = SegmentList.read('/home/albert.wandui/detchar'+\
    '/ER7/jul13/%s_ER7_segments.txt' %ifo)

# Read in the BBH triggers
bbh_trigs = SnglInspiralTable.read(bbhfile)
# We only want the triggers in the given segments
bbh_trigs = bbh_trigs.vetoed(segments)
#bbh_trigs.sort(key=lambda x: x.end_time + x.end_time_ns * 1.0e-9)

print "Read in all the BBH triggers!!!\n"
print "Let's start working on the Omicron triggers...\n"
# ---------------------------------------------------------------------------- #

# Read in all the Omicron caches
# Also get an idea of the speed of the code when reading from cache file vs
# letting vet get the data itself

Nchannels = len(channels)

def get_omicron_triggers(channel, ifo, segments, cachefile):
  print "Reading channel: %s\n" %channel
  with open(cachefile, 'r') as f:
    mycache = Cache.fromfile(f)
  # Let's try and catch failed reads
  try:
    triggers = get_triggers(ifo + ':' + channel, 'sngl_burst', segments,\
        cache=mycache)
  except:
    print "No Omicron triggers read for channel %s" %channel
    return None

  return triggers

# Use map to get all the data from the channels of interest at once
# Make lists of all the parameters
ifos = [ifo]*Nchannels
segs = [segments]*Nchannels

# Get all the cache files from the directory
cachefiles = glob.glob(os.path.join(omiccachedir, ifo + '/', '*.cache'))

# Now we can read in the triggers and also worry about performance
import time
t0 = time.time()
omic_trigger_tables = map(get_omicron_triggers, channels, ifos, segs, cachefiles)
t1  = time.time()
print "All the Omicron triggers for %d channels loaded." %len(omic_trigger_tables)
print "This took %f seconds to run to completion\n" %(t1 - t0)


window = 2.0
import h5py
f = h5py.File('vetosegments.hdf5', 'r')

snr = np.array(bbh_trigs.getColumnByName('snr')[:])
# ---------------------------------------------------------------------------- #
# First we will obtain the relevant metrics for all the channels and store them
# in a numpy array

# Metrics being considered
eff = get_metric('efficiency')
dt = get_metric('deadtime')
eff_over_dt = get_metric('efficiency/deadtime')
usep = get_metric('use percentage')
loudbysnr = get_metric('loudest event by snr')

mydtypes = [('channel', 'a50'), ('efficiency', float), ('deadtime', float),\
    ('efficiency_over_deadtime',float), ('use_percentage',float),\
    ('loudest_event', float)]

statistics = np.zeros((Nchannels), dtype=mydtypes)
# Make some basic histograms
for i in xrange(Nchannels):
  key = channels[i] +'/vetosegs'
  omic_trigs= omic_trigger_tables[i]
  vetosegs= SegmentList.read(f, key)
  after_trigs = bbh_trigs.veto(vetosegs)
  after_snr = np.array(after_trigs.getColumnByName('snr')[:])
  myflag = DataQualityFlag()
  myflag.active = vetosegs
  myflag.known = segments
  statistics[i] = (channels[i], eff(myflag, bbh_trigs).value, dt(myflag).value,\
      eff_over_dt(myflag, bbh_trigs).value, usep(myflag, omic_trigs).value,\
      loudbysnr(myflag, bbh_trigs).value)

  # Now do the plotting
  plt.figure(figsize=(12,10))
  bins = np.logspace(np.log10(5.5), np.log10(np.max(snr)), 50)
  labels = ['Before \nNTrigs=%d'%len(snr), 'After \nNTrigs=%d'%len(after_snr)]
  plt.hist([snr,after_snr], label=labels, bins=bins, histtype='step', log=True,\
      color=['b','r'])
  ax = plt.gca()
  ax.set_xlabel('Signal-to-noise ratio (SNR)')
  ax.set_ylabel('Counts (N)')
  ax.set_title('L1 Veto histogram for channel %s' %channels[i])
  plt.legend()
  plt.xlim(1,1e5)
  plt.grid(True, which="both")
  plt.savefig('%s_veto.png' %channels[i])
  plt.close()


## Write this data to a file
fmt = "%s %10.4f %10.4f %10.4f %10.4f %10.4f"
np.savetxt('vetostats.txt', statistics, fmt=fmt, delimiter=' ',\
    header='channel efficiency deadtime efficiency/deadtime use_percentage loudest_event')

print "All done!!!!"
f.close()
