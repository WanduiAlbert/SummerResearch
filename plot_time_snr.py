####################################
#import required 

from __future__ import (division, print_function)

import sys

from glue.segmentsUtils import fromsegwizard
from glue.segments import segment, segmentlist
from glue.lal import Cache
from glue.ligolw import table, utils, lsctables

from gwpy.detector import Channel
from gwpy.segments import Segment
from gwpy.table.lsctables import (SimBurstTable, SnglBurstTable)
from gwpy.plotter import (rcParams, HistogramPlot, EventTablePlot, SpectrumPlot)

from gwpy.segments import SegmentList

from pylal import ligolw_bucluster,snglcluster

from gwpy.toolkits.vet import get_triggers

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#Set the GPS start and end times
GPSstart =1117613088 
GPSend = 1117999425

### Parse the given arguments
# The argument is a path to the directory holding the cache files 
#to be processed

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('cachedir', action='store', type=str)
args = parser.parse_args()

#Use glob to generate a list of all the cache files
import glob
import os

#Make a list of all the cache files in the directory
cachelist = glob.glob(args.cachedir + '*.cache')

tag = ''

########## apply segment list (of good times)

segments = SegmentList.read('L1_ER7_segments.txt')

for cachefile in cachelist:
  ### open trigger cache
  # Make a tag for the saving the plot as well as the title
  # The tag is the name of the channel extracted from the path
  tag = cachefile.split('/')[-1]
  tag = tag.split('.')[0]
  tag = tag.replace('_Omicron','')
 
  print ('\n\nReading file: %s now ...\n' % tag)
  with open(cachefile, 'r') as fobj:
      cache = Cache.fromfile(fobj)


  ### read triggers
  # filter to select for triggers with frequency < 100
  #trigs = SnglBurstTable.read(cache, verbose=True, filt=lambda t: t.peak_frequency < 100)

  #filter to select for triggers with frequency <100 and snr <100
  trigs = get_triggers('L1:'+tag, 'sngl_burst', segments, cache=cache)


  ### check triggers read successfully
  if not trigs:
    print("    WARNING: No triggers for channel '%s'." % channel,
                file=sys.stderr)
  else:
      print('    %d triggers read' % len(trigs))

  ######## plot triggers and save plot

  plot = trigs.plot('time', 'snr', edgecolor='none')
  plot.set_ylabel('Central Frequency [Hz]')
  plot.set_yscale('log')
  title = r'L1:'+tag+' triggers'
  # Need to format _ to \_ for latex compatibility
  title  = title.replace('_', '{\_}')
  plot.set_title(title)
  # plot.add_colorbar(label='Signal-to-noise ratio')
  plt.savefig('%s_aux_snr_time.png' %tag)
