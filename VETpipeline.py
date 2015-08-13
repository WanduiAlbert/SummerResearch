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

# For each channel this defines the SNR thresholds to check
# for the 95,96,97,98,99,99.5% levels
thresholds = {"ASC-AS_B_RF36_I_YAW_OUT_DQ":[22.8,24.6,27.2,33.3,46.3,95.2],
"PEM-EY_MAG_EBAY_SUSRACK_QUAD_SUM_DQ":[18.4, 18.9, 19.5, 31.6, 116.9, 128.8],
"LSC-POP_A_LF_OUT_DQ":[43.6, 46.2, 49.9, 55.1,63.3,71.0],
"ASC-AS_B_RF45_Q_YAW_OUT_DQ":[468.6, 680.4, 1008.6, 1314.7, 7150.5, 8713.4],
"OMC-ASC_POS_Y_OUT_DQ":[16.8,17.6,18.6,19.9,22.0,23.8],
"PEM-EX_TILT_VEA_FLOOR_T_DQ":[19.6,20.3,21.0,22.4,25.5,29.4],
"PEM-EY_TILT_VEA_FLOOR_T_DQ":[19.6,20.1,20.9,22.0,24.5,27.9],
"OMC-PZT2_MON_AC_OUT_DQ":[729.2,10178,10767,10938,11105,11855],
"LSC-SRCL_OUT_DQ":[19.4, 20.4,21.8,24.4,31.1,44.3],
"IMC-REFL_DC_OUT_DQ":[46.5,50.8,57.0,67.6,95.8,130.7],
"ALS-X_REFL_ERR_OUT_DQ":[10.0,10.1,10.3,10.6,11.1,11.6],
"PEM-CS_MAG_EBAY_SUSRACK_Y_DQ":[25.8,26.0,26.3,26.6,27.6,35.3],
"ASC-Y_TR_B_NSUM_OUT_DQ":[20.8,21.2,21.5,21.8,25.3,51.2],
"ASC-AS_B_RF45_Q_PIT_OUT_DQ":[88.0, 315.8,570.5,1009.0,1951.8,8021.4],
"SUS-OMC_M1_ISIWIT_T_DQ":[18.4,19.6,23.7,27.7,40.3,426.9],
"PSL-ISS_AOM_DRIVER_MON_OUT_DQ":[16.6,24.1,28.2,20.8,34.2,38.3],
"LSC-PRCL_OUT_DQ":[21.3,22.6,24.3,26.8,30.2,35.6]}

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
peaktime_all_channels = map(lambda x: \
    np.array(x.getColumnByName('peak_time')[:], dtype=float) + \
    np.array(x.getColumnByName('peak_time_ns')[:], dtype=float) * 1.0e-9,\
    omic_trigger_tables)
#snr = np.array(map(lambda x: x.getColumnByName('snr')[:], omic_trigger_tables)

# Calculate the offsets for a single channel. omic_times is the array of all the
# peaktimes for a single channel and bbhtimes is the array of endtimes for the
# BBH triggers. At the end we take the transpose so as to ensure that the shape
# of the resulting array is (NumOmicron, NumBBH). Each column represents the
# offsets between the Omicron triggers and a single BBH trigger end time.
window = 2.0
import h5py
f = h5py.File('vetosegments.hdf5', 'w')

# Function that gets the veto time for a single bbhtime
def get_vetotimes(omic_peaktimes, end_times, veto_segs):
  for endtime in end_times:
    if np.any(np.abs(endtime - omic_peaktimes) <= window/2.0):
      veto_segs.append((float(endtime - window/2.0), float(endtime + window/2.0)))

percentile = [95,96,97,98,99,99.5]

for j in xrange(6):
  print "Get the offsets for the %d % percentile.\n" %percentile[j]
  thresh = f.create_group('%.1f' %percentile[j])
  for i in xrange(len(channels)):
    print "Working on channel %d now...\n" %i
    # Get the SNR threshold to use for this channel.
    snr_thresh = thresholds[channels[i]][j]
    omic_peaktimes = peaktime_all_channels[i] # one channel at a time!
    selection = omic_peaktimes.copy()
    selection.extend(filter(lambda x: x.snr > snr_thresh,omic_peaktimes))
    veto_segs = []
    t0 = time.time()
    get_vetotimes(selection, end_times, veto_segs)
    t1 = time.time()
    print "This took %f seconds to run to completion\n" %(t1 - t0)

    print "Offsets calculated. Coalesce the segments now\n"
    veto_segs = SegmentList(veto_segs)
    # Merge contiguous veto sections and sort the list of segments
    veto_segs.coalesce()

    # Write all the segments to disk. Read in the [0, 2] columns to recover the data
    # We will use h5py to write out all the segments in groups organized by the name
    # of the channel
    print "Write the segments to file\n"
    grp = thresh.create_group('%s' %channels[i])
    SegmentList.write(veto_segs, grp, 'vetosegs')
    print "All done!"

f.close()
print "Finished computing the offsets for all the channels!!!! Wooohoo!!!\n"
