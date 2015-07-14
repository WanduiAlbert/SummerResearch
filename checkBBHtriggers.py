#! /usr/bin/env python

# This is a simple script to check if all the data in the BBH trigger
# files is according to specification. Particularly, that all the triggers 
# are from the same detector and that no triggers overlap within a 0.1 second
# window

import sys
import numpy as np
from gwpy.table.lsctables import SnglInspiralTable

triggers = sys.argv[1]
gwdata = SnglInspiralTable.read(triggers)

ifos = np.array(gwdata.getColumnByName('ifo')[:])
same_ifo = np.any(ifos != ifos[0])

if same_ifo:
  sys.exit("Error: Not all triggers are from the same detector!")

etimes = np.array(gwdata.getColumnByName('end_time')[:])
etimes_ns = np.array(gwdata.getColumnByName('end_time_ns')[:])

etimes = etimes + etimes_ns * 1.0e-9

etimes = np.sort(etimes)
delta_t = np.diff(etimes)

overlap = np.any(delta_t < 0.1)

if overlap:
  sys.exit("Error: Some triggers overlap to within 0.1 seconds!")

print "BBH trigger file is according to specification!"
