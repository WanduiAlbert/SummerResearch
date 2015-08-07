# Import statements

from __future__ import division

import sys

import numpy as np

from gwpy.segments import Segment, SegmentList
from gwpy.table.lsctables import SnglBurstTable, SnglInspiralTable

from glue.lal import Cache

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
omic_triggers = map(get_omicron_triggers, channels, ifos, segs, cachefiles)
t1  = time.time()
print "All the Omicron triggers for %d channels loaded." %len(omic_triggers)
print "This took %f seconds to run to completion\n" %(t1 - t0)

# Get the peak times of the Omicron triggers
omic_peaktimes = np.array(map(lambda x: x.get_peak(), omic_triggers))

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
tiledomictimes = np.array(map(lambda x: np.tile(x, NumBBH), omic_peaktimes))
tiledbbhtimes = end_times.reshape(NumBBH, 1)

# Get all the offsets for all the channels at one time.
offsets = np.array(map(lambda x: np.abs(tiledbbhtimes - x), tiledomictimes))
print "This worked out fine thus far!!!"



