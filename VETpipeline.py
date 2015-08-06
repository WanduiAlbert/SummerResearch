# Import statements

from __future__ import division

import sys

import numpy as np

from gwpy.segments import Segment, SegmentList
from gwpy.table.lsctables import SnglBurstTable, SnglInspiralTable

from glue.lal import Cache

from gwpy.toolkits.vet import (get_triggers,get_segments,get_metric)

# Save all the channels in a python script
from channels import channels

ifo = sys.argv[1]
bbhfile = sys.argv[2]
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
end_times = bbh_trigs.get_end()

m1 = bbh_trigs.getColumnByName('mass1')
m2 = bbh_trigs.getColumnByName('mass2')

M = m1+m2
eta = (m1*m2)/M**2
mchirp = eta**(3./5)*M

# Consider appending the Mchirp column to the SnglInspiralTable
# ---------------------------------------------------------------------------- #

# Read in all the Omicron caches
# Also get an idea of the speed of the code when reading from cache file vs
# letting vet get the data itself

Nchannels = len(channels)

def get_omicron_triggers(channel, ifo, segments, cachefile):
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
cachefiles = [cachefile]*Nchannels

import time
t0 = time.time()
omic_triggers = map(get_omicron_triggers, channels, ifos, segs, cachefiles)
t1  = time.time()
print "All the Omicron triggers for %d channels loaded." %len(omic_triggers)
print "This took %f seconds to run to completion" %(t1 - t0)

# Get the peak times of the Omicron triggers
omic_endtimes = map(lambda x: x.get_peak(), omic_triggers)

# ---------------------------------------------------------------------------- #
# Now let's calculate the offsets between the BBH and Omicron triggers 

