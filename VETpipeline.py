# Import statements

import sys

import numpy as np

from gwpy.segments import Segment, SegmentList
from gwpy.table.lsctables import SnglBurstTable, SnglInspiralTable

from glue.lal import Cache

from  gwpy.time import LIGOTimeGPS

from gwpy.toolkits.vet import (get_triggers,get_segments,get_metric)

import glob, os

# Save all the channels in a python script
channels = ["ALS-X_REFL_ERR_OUT_DQ",
"ASC-AS_B_RF36_I_YAW_OUT_DQ",
"ASC-AS_B_RF45_Q_PIT_OUT_DQ",
"ASC-AS_B_RF45_Q_YAW_OUT_DQ",
"ASC-Y_TR_B_NSUM_OUT_DQ",
"IMC-REFL_DC_OUT_DQ",
"PSL-ISS_AOM_DRIVER_MON_OUT_DQ"]

#channels = ["ASC-AS_B_RF36_I_YAW_OUT_DQ",
#"PEM-EY_MAG_EBAY_SUSRACK_QUAD_SUM_DQ",
#"LSC-POP_A_LF_OUT_DQ",
#"ASC-AS_B_RF45_Q_YAW_OUT_DQ",
#"OMC-ASC_POS_Y_OUT_DQ",
#"PEM-EX_TILT_VEA_FLOOR_T_DQ",
#"PEM-EY_TILT_VEA_FLOOR_T_DQ",
#"OMC-PZT2_MON_AC_OUT_DQ",
#"LSC-SRCL_OUT_DQ",
#"IMC-REFL_DC_OUT_DQ",
#"ALS-X_REFL_ERR_OUT_DQ",
#"PEM-CS_MAG_EBAY_SUSRACK_Y_DQ",
#"ASC-Y_TR_B_NSUM_OUT_DQ",
#"ASC-AS_B_RF45_Q_PIT_OUT_DQ",
#"SUS-OMC_M1_ISIWIT_T_DQ",
#"PSL-ISS_AOM_DRIVER_MON_OUT_DQ",
#"LSC-PRCL_OUT_DQ"]

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

# We need to extract the chirp mass and the end times for these triggers
end_times = np.array(bbh_trigs.getColumnByName('end_time')[:], dtype=float) +\
    np.array(bbh_trigs.getColumnByName('end_time_ns')[:], dtype=float) * 1.0e-9

m1 = np.array(bbh_trigs.getColumnByName('mass1')[:])
m2= np.array(bbh_trigs.getColumnByName('mass2')[:])

M = m1+m2
eta = (m1*m2)/M**2
mchirp = eta**(3./5)*M

print "Read in all the BBH triggers!!!\n"
print "Let's start working on the Omicron triggers...\n"
# Consider appending the Mchirp column to the SnglInspiralTable
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

# ---------------------------------------------------------------------------- #
# Now let's calculate the offsets between the BBH and Omicron triggers.
# For each channel, we will make a two-dimensional array of repeated Omicron
# endtimes. Then we can use vector operations and completely avoid loops in this
# code. The BBH_endtimes have to be transposed in order for this to work. We
# will then check to see if there are any Omicron triggers in the coincidence
# window of a BBH trigger and get the segment times for those BBH triggers so as
# to veto them.

# Simply preparing the data for vectorized computing of the offsets
print "Starting the computation of the offsets...\n"
NumBBH = len(end_times)

# Here is where I should define the threshold on the SNR
# Get the peak times of the Omicron triggers

#snr = np.array(map(lambda x: x.getColumnByName('snr')[:], omic_trigger_tables)

# Calculate the offsets for a single channel. omic_times is the array of all the
# peaktimes for a single channel and bbhtimes is the array of endtimes for the
# BBH triggers. At the end we take the transpose so as to ensure that the shape
# of the resulting array is (NumOmicron, NumBBH). Each column represents the
# offsets between the Omicron triggers and a single BBH trigger end time.
window = 2.0
import h5py
f = h5py.File('vetosegments.hdf5', 'w')

def filter_threshold(omic_trigs,snr_thresh):
  selection = omic_trigs.copy()
  selection.extend(filter(lambda x: x.snr > snr_thresh,omic_trigs))
  return selection

def get_percentile(before, level):
  """
  Get a percentile level for the snr of the triggers from a channel.
  This is useful in evaluating appropriate levels at which to set the
  snr threshold.
  """
  snr = np.array(before.getColumnByName('snr')[:])
  return np.percentile(snr, level)

# Function that gets the veto time for a single bbhtime
def get_vetotimes(omic_peaktimes, end_times, veto_segs):
  for endtime in end_times:
    if np.any(np.abs(endtime - omic_peaktimes) <= window/2.0):
      veto_segs.append((float(endtime - window/2.0), float(endtime + window/2.0)))

# Percentiles at which to check our data
percentile = np.arange(0,100,2)

# Loop over the percentiles and extract the relevant statistics
for j in xrange(len(percentile)):
  print "Get the offsets for the %d %% percentile.\n" %percentile[j]
  thresh_grp = f.create_group('%.1f' %percentile[j])

  for i in xrange(len(channels)):
    print "Working on channel %d now...\n" %i
    # Get the SNR threshold to use for this channel.
    omic_trigs = omic_trigger_tables[i] # get triggers from one channel at a time!
    snr_thresh = get_percentile(omic_trigs, percentile[j])
    # Apply the snr filter to the triggers and get their peaktimes
    selection = filter_threshold(omic_trigs, snr_thresh)
    peaktime = np.array(selection.getColumnByName('peak_time')[:], dtype=float) + \
      np.array(selection.getColumnByName('peak_time_ns')[:], dtype=float) * 1.0e-9
    # Now we can calculate the offsets for this channel
    veto_segs = []
    t0 = time.time()
    get_vetotimes(peaktime, end_times, veto_segs)
    t1 = time.time()
    print "This took %f seconds to run to completion\n" %(t1 - t0)

    # Coalesce the segments and write them out to the file
    print "Offsets calculated. Coalesce the segments now\n"
    veto_segs = SegmentList(veto_segs)
    # Merge contiguous veto sections and sort the list of segments
    veto_segs.coalesce()

    # Write all the segments to disk. Read in the [0, 2] columns to recover the data
    # We will use h5py to write out all the segments in groups organized by the name
    # of the channel
    print "Write the segments to file\n"
    grp = thresh_grp.create_group('%s' %channels[i])
    SegmentList.write(veto_segs, grp, 'vetosegs')
    print "All done!"

f.close()
print "Finished computing the offsets for all the channels!!!! Wooohoo!!!\n"
