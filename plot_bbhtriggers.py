#! /usr/bin/env python

from __future__ import (division, print_function)

import sys

from gwpy.segments import SegmentList
from gwpy.table.lsctables import SnglInspiralTable

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

bbh_file = sys.argv[1]

trigs = SnglInspiralTable.read(bbh_file)
segs = SegmentList.read('L1_ER7_segments.txt')
trigs = trigs.vetoed(segs)

plot = trigs.plot('time', 'snr', edgecolor='none')#, epoch=1117378816)
#plot.set_xlim(1117378816, 1117378816+(24*3600*11.0))
plot.set_ylabel('SNR')
plot.set_yscale('log', nonposy='clip')
plot.set_title('BBH triggers during the ER7 run')
plot.savefig('H1_BBH_SNR.png')


