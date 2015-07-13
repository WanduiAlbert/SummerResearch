#! /bin/bash

SITE='L1' #change this to H1 if working in the Hanford clusteri
CHANNEL='GDS-CALIB_STRAIN_Omicron' # GSD calib strain

ls /home/detchar/triggers/ER7/$SITE/$CHANNEL/*/* | lalapps_path2cache > /home/albert.wandui/caches/hoft_Omicron_ER7.cache 
