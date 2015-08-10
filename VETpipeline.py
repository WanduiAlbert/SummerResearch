# Import statements

from __future__ import division

import sys

import numpy as np

from gwpy.segments import Segment, SegmentList
from gwpy.table.lsctables import SnglBurstTable, SnglInspiralTable

from glue.lal import Cache

from  gwpy.time import LIGOTimeGPS

from gwpy.toolkits.vet import (get_triggers,get_segments,get_metric)

import glob, os
# Save all the channels in a python script
from channels import channels

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
end_times = np.array(bbh_trigs.get_end())

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
  with open(cachefile) as f:
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
peaktime_all_channels = map(lambda x: x.get_peak(), omic_trigger_tables)
#snr = np.array(map(lambda x: x.getColumnByName('snr')[:], omic_trigger_tables)

# Calculate the offsets for a single channel. omic_times is the array of all the
# peaktimes for a single channel and bbhtimes is the array of endtimes for the
# BBH triggers. At the end we take the transpose so as to ensure that the shape
# of the resulting array is (NumOmicron, NumBBH). Each column represents the
# offsets between the Omicron triggers and a single BBH trigger end time.
#def get_offset(omic_times, end_times):
#  print "I am currently working on channel %d" %pos
#  return map(lambda omictime: np.abs(end_times - omictime),\
#      omic_times)
#
#print "Now lets tile the Omicron end times and calculate the offsets.\n"
#
#for i in xrange(len(peaktime_all_channels[0])):
#  print "Index : %d\n" %i
#  offsets = get_offset(peaktime_all_channels[0], end_times)
#
#print offsets.shape
##allofssets = map(get_offset, peaktime_all_channels[0],\ [end_times])#*Nchannels)
#
#print "This worked out fine thus far!!!"
window = 2.0
omic_peaktimes = peaktime_all_channels[0] # one channel at a time!
# Function that gets the veto time for a single bbhtime
def get_vetotimes(bbhtime, omic_peaktimes):
  vetosegs = []
  # It is interesting to note that omic_peaktimes - bbhtime works but
  # bbhtime-omic_peaktimes doesn't. This symmetry is broken by the fact that
  # the subtraction operator for LIGOTimeGPS is already overloaded. So I chose
  # to do the omic-bbh since we will be taking the absolute value.
  # Alternatively, you could np.array(bbhtime).
  if np.any(np.abs(omic_peaktimes - bbhtime) <= window/2.0):
    vetosegs += [Segment(bbhtime-window/2.0, bbhtime+window/2.0)]
  return vetosegs

all_veto_segs = []
for bbhtime in end_times:
  all_veto_segs += get_vetotimes(bbhtime, omic_peaktimes)
all_veto_segs = SegmentList(all_veto_segs)
# Merge contiguous veto sections and sort the list of segments
all_veto_segs.coalesce()

# Write all the segments to disk. Read in the [0, 2] columns to recover the data
# We will use h5py to write out all the segments in groups organized by the name
# of the channel
import h5py
f = h5py.File('vetosegments.hdf5', 'w')
grp = f.create_group('%s' %channels[0])
SegmentList.write(all_veto_segs, grp, 'vetosegs')
f.close()
