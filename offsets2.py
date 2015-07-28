#! /usr/bin/env python

# This is a script to generate time offsets between BBH triggers and Omicron
# triggers. A histogram of the time offsets is plotted as a function of the 
# chirp mass of the BBH trigger.

############################################################################
#Imports

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from gwpy.table.lsctables import (SnglInspiralTable, SnglBurstTable)
from gwpy.segments import SegmentList
from glue.lal import Cache

import sys
import os

###############################################################################
# This script has been split into functions since the read omic trigger function
# takes a really long time but only needs to be called once. Best run from
# ipython

# Load the BBH triggers into a SnglInspiralTable
def load_bbh_trigs(bbhtrigfile, segs):
  # Read in the BBh triggers
  bbh_trigs = SnglInspiralTable.read(bbhtrigfile)

  # Check if BBH triggers have been read in successfully
  if not bbh_trigs:
    sys.exit("ERROR: No triggers for BBH file: %s" % bbhtrigfile)
  else:
    print "%d BBH triggers read" % len(bbh_trigs)

  #Get the BBH triggers that lie within the valid segment list
  bbh_trigs = bbh_trigs.vetoed(segs)
  # Sort the bbh_triggers by their end times. This will be useful later
  bbh_trigs.sort(key=lambda x: x.end_time + x.end_time_ns * 1.0e-9)
  return bbh_trigs

# Load the Omicron triggers from the cache file into a SnglBurstTable
# This call takes a long time to execute
def load_omic_trigs(omicroncachefile, segs):
  # Read in the Omicron triggers
  with open(omicroncachefile, 'r') as cachefile:
    cache = Cache.fromfile(cachefile)

  omic_trigs = SnglBurstTable.read(cache, verbose=True, filt=lambda x: x.snr <\
      100 and x.peak_frequency < 100)

  # Check if Omicron triggers have been read in successfully
  if not omic_trigs:
    sys.exit("ERROR: No triggers for Omicron channel: %s" % cachefile.split('.')[0])
  else:
    print "%d Omicron triggers read" % len(omic_trigs)

  #Get the Omicron triggers that lie within the valid segment list
  omic_trigs = omic_trigs.vetoed(segs)
  return omic_trigs

# Extract the end times and calculate the chirp mass for the BBH triggers
def get_bbh_params(bbh_trigs):
  bbh_endtimes = np.array(bbh_trigs.getColumnByName('end_time')[:]) +\
      np.array(bbh_trigs.getColumnByName('end_time_ns')[:]) * 1.0e-9
  mass1 = np.array(bbh_trigs.getColumnByName('mass1')[:])
  mass2 = np.array(bbh_trigs.getColumnByName('mass2')[:])
  M = mass1+ mass2
  eta = mass1*mass2/M**2
  mchirp = eta**(3./5)*M
  return bbh_endtimes, mchirp

# Extract the peak times for the Omicron triggers
def get_omic_params(omic_trigs):
  omic_times = np.array(omic_trigs.getColumnByName('peak_time')[:]) + \
      np.array(omic_trigs.getColumnByName('peak_time_ns')[:]) * 1.0e-9
  return omic_times

# For each BBH trigger, gets the Omicron triiggers within a window around 
# the end time of the trigger and returns an array of the offsets of the
# triggers.
def get_offsets(bbh_endtimes, mchirp,  omic_times, window):
  # Step through the BBH triggers and check for the omicron triggers in the
  # coincidence window.
  all_trigs = {}
  for time, mch in zip(bbh_endtimes, mchirp):
    diff = omic_times - time
    # In order to check for overlaps, need to get the indices at which offsets 
    # takes values. This indices will be used to make a dict that keeps track of
    # the omicron triggers already in offset. If an omicron trigger already in
    # this dict is closer to another BBH trigger than the current one, then its
    # offsets value is replaced, otherwise it is not.
    offset_index = np.where(np.abs(diff) <= (window/2.0))[0]
    # Hopefully, we can loop over these indices since they should be much fewer
    # in number than the full set of omic_triggers. Also store the chirp mass
    # associated with the bbh time.
    for index in offset_index:
      # add the index to the dict since we haven't seen it yet.
      if index not in all_trigs:
        all_trigs[index] = (diff[index], mch)
      else:
        # if we have already seen the index then ensure that we associate it
        # with the closest possible trigger
        if np.abs(all_trigs[index][0]) > np.abs(diff[index]):
          all_trigs[index] = (diff[index], mch)


  all_offsets = np.array(all_trigs.values(), dtype=[('offset',\
      np.float32),('mchirp', np.float32)])
  all_offsets.sort(order='mchirp')
  return all_offsets

# Simple script to load the segment lists for valid data from a text file
def load_segs():
  segs = SegmentList.read('L1_ER7_segments.txt')
  return segs

# Simple script to generate the tag for legends and plot names
def get_tag(omicroncachefile):
  tag = omicroncachefile.split('/')[-1]
  tag = tag.split('.')[0]
  return tag

# This makes subplots of the offsets in groupings of the chirp
# mass. The hope is that this will make the trends more clear
# and the patterns way more distinguishable
def plot_subplots(offsets, tag, n, window):
  fig, axarr = plt.subplots(4, sharex=True)
  sel1 = np.logical_and(offsets['mchirp'] > 0.0, offsets['mchirp'] <= 5.0)
  sel2 = np.logical_and(offsets['mchirp'] > 5.0, offsets['mchirp'] <= 10.0)
  sel3 = np.logical_and(offsets['mchirp'] > 10.0, offsets['mchirp'] <= 15.0)
  sel4 = np.logical_and(offsets['mchirp'] > 15.0, offsets['mchirp'] <= 100.0)
  
  sel = [sel1, sel2, sel3, sel4]
  labels = [r'$\mathcal{M}\ \leq\ 5M_{\odot}$', r'$5M_{\odot}\ <\ \mathcal{M}\ \leq\ 10M_{\odot}$'\
      ,r'$10M_{\odot}\ <\ \mathcal{M}\ \leq\ 15M_{\odot}$', r'$\mathcal{M}\ >\ 15M_{\odot}$']
  
  for i in xrange(len(axarr)):
    if not np.any(sel[i]):
      continue

    axarr[i].hist(offsets['offset'][sel[i]], histtype='step', bins=n, label=labels[i])
    axarr[i].set_yscale('log', nonposy='clip')
    axarr[i].set_xlim(-0.25, 0.25)
    axarr[i].set_ylabel(r'N')
    axarr[i].grid(True)
    axarr[i].legend(loc="upper left", bbox_to_anchor=(1,1))

  axarr[3].set_xlabel(r'Offset  [Sec]')
  fig.subplots_adjust(hspace=0.5)
  plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

  axarr[0].set_title( tag)
  plt.savefig(tag + '_subplotshistogram.png', bbox_inches='tight')
  #, bbox_extra_artists=(lgd,), bbox_inches='tight')
  plt.close()

if __name__=='__main__':
  segs = load_segs()
  n = 50
  window = 0.5 
  tag = ''
  bbhfile = sys.argv[1]
  bbh_trigs = load_bbh_trigs(bbhfile, segs)
  bbh_endtimes, mchirp = get_bbh_params(bbh_trigs)
  
  omiccachedir = sys.argv[2]
  import glob
  omiccachelist = glob.glob(omiccachedir + '*.cache')
  plt.figure()
  
  for omiccachefile in omiccachelist:
    #omiccachefile = omiccachelist[i]
    print "\n\nCurrently reading file: %s \n" % omiccachefile
    omic_trigs = load_omic_trigs(omiccachefile, segs)
    tag = get_tag(omiccachefile)
    print "\n\nExtracting the peak times ...\n"
    omic_times = get_omic_params(omic_trigs)
    print "\nDone...\n Calculating the offsets...\n"
    offsets = get_offsets(bbh_endtimes, mchirp, omic_times, window)
    print "\n Done... \n Adding the mean offsets to the plot ...\n"
    print "\n\n Done. Moving on to the next file in the list. \n"
    plot_subplots(offsets, tag, n, window)

  print "\nAll plots are completed!\n"
