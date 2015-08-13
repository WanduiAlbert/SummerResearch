# Import statements
from __future__ import division

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

from gwpy.plotter import EventTablePlot

from gwpy.toolkits.vet import Metric, register_metric
from astropy.units import Unit
from astropy.io import ascii

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



# Redefine use percentage
def my_use_percentage(segments,before,after=None):
  """
  Calculate the use percentage as the ratio of the triggers in
  an auxiliary channel used to the total number of triggers in
  the channel. The standard definition checks segments instead.
  """
  return len(before.vetoed(segments.active))/len(before)
register_metric(Metric(my_use_percentage, "my use percentage", unit=Unit('%')))

# Defining the functions for all 9 plots and the summary statistics
def summary_stats(statistics, bbh_trigs, omicron_trigs, channel, vetosegs, segments):
  # Metrics being considered
  eff = get_metric('efficiency')
  dt = get_metric('deadtime')
  eff_over_dt = get_metric('efficiency/deadtime')
  usep = get_metric('my use percentage')
  loudbysnr = get_metric('loudest event by snr')
  myflag = DataQualityFlag()
  myflag.active = vetosegs
  myflag.known = segments
  statistics[i] = (channel, eff(myflag, bbh_trigs).value, dt(myflag).value,\
      eff_over_dt(myflag, bbh_trigs).value, usep(myflag, omic_trigs).value,\
      loudbysnr(myflag, bbh_trigs).value)

def downtime(vetosegs, segments, channel):
  myflag = DataQualityFlag()
  myflag.active = vetosegs
  myflag.known = segments
  plot = myflag.plot()
  plot.set_title(r'Active and vetoed segments for %s' %channel)
  plot.savefig(r'%s_downtime.png'%channel)

def histogram(bbh_trigs, vetosegs, channel):
  snr = np.array(bbh_trigs.getColumnByName('snr')[:])
  after_trigs = bbh_trigs.veto(vetosegs) # Triggers that remained
  after_snr = np.array(after_trigs.getColumnByName('snr')[:])
  plt.figure(figsize=(12,10))
  bins = np.logspace(np.log10(5.5), np.log10(np.max(snr)), 50)
  labels = [r'Before \nNTrigs=%d'%len(snr), r'After \nNTrigs=%d'%len(after_snr)]
  plt.hist([snr,after_snr], label=labels, bins=bins, histtype='step',\
      color=['b','r'])
  ax = plt.gca()
  ax.set_xlabel('Signal-to-noise ratio (SNR)')
  ax.set_ylabel('Counts (N)')
  ax.set_yscale('log', nonposy='clip')
  ax.set_xscale('log', nonposx='clip')
  ax.set_title(r'%s histogram for channel %s' %(ifo, channel))
  plt.legend()
  plt.xlim(5.5, np.max(snr))
  plt.grid(True, which="both")
  plt.savefig(r'%s_histogram.png' %channel)
  plt.close()

def time_snr(bbh_trigs, vetoed_trigs, channel):
  labels = [r"All %d" %len(bbh_trigs), r"Vetoed %d" %len(vetoed_trigs)]
  plot = bbh_trigs.plot('time', 'snr',\
      edgecolor='none', label=labels[0])
  plot.add_table(vetoed_trigs, 'time', 'snr',\
      edgecolor='none', label=labels[1],c='r')
  plot.set_ylabel('SNR')
  plot.set_yscale('log',nonposy='clip')
  plot.set_title(r'Detector=%s, Veto Channel=%s' %(ifo, channel))
  ax = plot.gca()
  lgd = ax.legend(loc="upper left", bbox_to_anchor=(1,1))
  save = channel.replace('{\_}', '_')
  plot.savefig(r'%s_main_time_snr.png' %save,\
    bbox_extra_artists=(lgd,), bbox_inches='tight')

def aux_time_freq(omic_trigs, vetoed_omic_trigs, channel):
  labels = [r"All %d" %len(omic_trigs), r"Used %d" %len(vetoed_omic_trigs)]
  plot = omic_trigs.plot('time', 'central_freq',\
      edgecolor='none', label=labels[0])
  plot.add_table(vetoed_omic_trigs, 'time', 'central_freq',\
      edgecolor='none', label=labels[1],c='r')
  plot.set_ylabel('Central Frequency [Hz]')
  plot.set_title(r'Detector=%s, Veto Channel=%s' %(ifo, channel))
  plot.set_yscale('log',nonposy='clip')
  ax = plot.gca()
  lgd = ax.legend(loc="upper left", bbox_to_anchor=(1,1))
  save = channel.replace('{\_}', '_')
  plot.savefig(r'%s_aux_time_freq.png' %save,\
    bbox_extra_artists=(lgd,), bbox_inches='tight')

def aux_snr_freq(omic_trigs, vetoed_omic_trigs, channel):
  labels = [r"All %d" %len(omic_trigs), r"Used %d" %len(vetoed_omic_trigs)]
  plot = omic_trigs.plot('central_freq','snr',\
      edgecolor='none', label=labels[0])
  plot.add_table(vetoed_omic_trigs,'central_freq','snr',\
      edgecolor='none', label=labels[1],c='r')
  plot.set_xlabel('Central Frequency [Hz]')
  plot.set_ylabel('SNR')
  plot.set_title(r'Detector=%s Veto Channel=%s' %(ifo, channel))
  plot.set_yscale('log',nonposy='clip')
  plot.set_yscale('log', nonposx='clip')
  ax = plot.gca()
  lgd = ax.legend(loc="upper left", bbox_to_anchor=(1,1))
  save = channel.replace('{\_}', '_')
  plot.savefig(r'%s_aux_freq_snr.png' %save,\
    bbox_extra_artists=(lgd,), bbox_inches='tight')

def aux_snr_time(omic_trigs,vetoed_omic_trigs, channel):
  labels = [r"All %d" %len(omic_trigs), r"Used %d" %len(vetoed_omic_trigs)]
  plot = EventTablePlot(omic_trigs, 'time', 'snr',\
      edgecolor='none', label=labels[0])
  plot.add_table(vetoed_omic_trigs, 'time', 'snr',\
      edgecolor='none', label=labels[1],c='r')
  plot.set_ylabel('SNR')
  plot.set_yscale('log',nonposy='clip')
  plot.set_title(r'Detector=%s, Veto Channel=%s' %(ifo, channel))
  ax = plot.gca()
  lgd = ax.legend(loc="upper left", bbox_to_anchor=(1,1))
  save = channel.replace('{\_}', '_')
  plot.savefig(r'%s_aux_time_snr.png' %save,\
    bbox_extra_artists=(lgd,), bbox_inches='tight')

def cumulative_histogram(omic_trigs,vetosegs,channel):
  plot = omic_trigs.hist(snr, log=True, logbins=True, histtype='stepfilled',\
    color='g', cumulative=True, normed=True)
  plot.set_xlabel('Signal-to-noise ratio (SNR)')
  plot.set_ylabel('Normed Counts')
  plot.set_yscale('log', nonposy='clip')
  plot.set_xscale('log', nonposx='clip')
  plot.set_title(r'%s Cumulative histogram for channel %s' %(ifo, channel))
  plot.savefig(r'%s_cumulative_histogram.png' %channel)

mydtypes = [('channel', 'a50'), ('efficiency', float), ('deadtime', float),\
      ('efficiency/deadtime',float), ('use_percentage',float),\
      ('loudest_event', float)]
statistics = np.zeros((Nchannels), dtype=mydtypes)
# Make some basic histograms
for i in xrange(Nchannels):
  key = channels[i] +'/vetosegs'
  omic_trigs= omic_trigger_tables[i]
  vetosegs= SegmentList.read(f, key)
  vetoed_trigs= bbh_trigs.vetoed(vetosegs)
  vetoed_omic_trigs= omic_trigs.vetoed(vetosegs) #Triggers that were vetoed
  summary_stats(statistics, bbh_trigs, omic_trigs, channels[i], vetosegs, segments)
  cumulative_histogram(omic_trigs,vetosegs,channels[i])
  # Now do the plotting
  # histogram(bbh_trigs, vetosegs, channels[i])
  # channel= channels[i]
  # channel= channel.replace('_','{\_}')
  # print "Working on the time snr plot now \n"
  # time_snr(bbh_trigs, vetoed_trigs, channel)
  # print "Working on the time frequency plot for the aux channel now \n"
  # aux_time_freq(omic_trigs, vetoed_omic_trigs, channel)
  # print "Working on the snr frequency plot for the aux channel now \n"
  # aux_snr_freq(omic_trigs, vetoed_omic_trigs, channel)
  # print "Working on the snr time plot for the aux channel now \n"
  # aux_snr_time(omic_trigs,vetoed_omic_trigs, channel)
  # downtime(vetosegs, segments, channel)
  print "Done! Moving on to the next channel!!!!\n"

## Write this data to a file
# We first sort our data by efficiency/deadtime, then efficiency
# and finally deadtime
statistics[::-1].sort(order=["efficiency/deadtime","deadtime","efficiency"])
fmt = {"channel":"%-50s","efficiency":"%10.4f", "deadtime":"%10.4f",\
  "efficiency/deadtime":"%10.4f","use percentage":"%10.4f",'loudest event by snr':"%10.4f"}
names = ["channel", "efficiency", "deadtime",\
  "efficiency/deadtime", "use percentage","loudest event by snr"]
ascii.write(statistics, output="vetostats.txt",format="fixed_width", names=names,\
  comment="#", formats=fmt, delimiter="|", delimiter_pad=" ")

print "All done!!!!"
f.close()
