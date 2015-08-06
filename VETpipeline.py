# Import statements

from __future__ import division

import sys

import numpy as np

from gwpy.segments import Segment, SegmentList
from gwpy.table.lsctables import SnglBurstTable, SnglInspiralTable

from glue.lal import Cache

from gwpy.toolkit.vet import (get_triggers,get_segments,get_metric)

ifo = sys.argv[1]
bbhfile = sys.argv[2]
omiccachedir = sys.argv[3]

# Read in the segment file
segments = SegmentList.read('/home/albert.wandui/detchar/ER7/jul13/%s_ER7_segments.txt' %ifo)

# Read in the BBH triggers
bbh_trigs = SnglInspiralTable.read(bbhfile)
# We only want the triggers in the given segments
bbh_trigs = bbh_trigs.vetoed(segments)

# We need to extract the chirp mass and the end times for these triggers
end_times = bbh_trigs.get_end()

m1 = bbh_trigs.getColumnByName('mass1')
m2 = bbh_trigs.getColumnByName('mass2')
mchirp = (m1*m2)/(m1+m2)**2


