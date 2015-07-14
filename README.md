# SummerResearch

This is the repository that holds all the code that I will write/use that
pertains to the analysis of BBH triggers and Omicron triggers for aLIGO during
the ER7 run.

Current scripts and their descriptions:
  1. plot_alltriggers: reads in BBH triggers in xml or npz format and generates
     histograms of the SNR. The SNR histograms, logarithmic in the number of 
     events, show a general straight line trend, characteristic of a gaussian 
     distribution with a non-gaussian tail at high SNR due to glitches. That
     this is so is not very clear with the limited data from ER7.
  
  2. makecache.sh: Generates cache files for a list of relevant channels. For L1
     the channels chosen are the top contenders as determined using Hveto here:
     https://ldas-jobs.ligo-la.caltech.edu/~areeda/L1-omicron_BOTH-1117400416-842400-DARM/
     The cache files are then read into plot_time_snr and other related scripts.

  3. plot_time_snr.py:Takes in as an argument the path to the folder containing
     a bunch of cache files from the Omicron triggers. These files are read and
     the data from the triggers extracted. So far, plots of the time vs central 
     frequency and SNR have been made. Other relevant plots will be made as need
     be.


