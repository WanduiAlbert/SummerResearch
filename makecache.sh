#! /bin/bash

SITE='H1' #change this to H1 if working in the Hanford clusteri

CHANNELS="ASC-AS_B_RF36_I_YAW_OUT_DQ_Omicron
PEM-EY_MAG_EBAY_SUSRACK_QUAD_SUM_DQ_Omicron
LSC-POP_A_LF_OUT_DQ_Omicron
ASC-AS_B_RF45_Q_YAW_OUT_DQ_Omicron
OMC-ASC_POS_Y_OUT_DQ_Omicron
PEM-EX_TILT_VEA_FLOOR_T_DQ_Omicron
PEM-EY_TILT_VEA_FLOOR_T_DQ_Omicron
OMC-PZT2_MON_AC_OUT_DQ_Omicron
LSC-SRCL_OUT_DQ_Omicron
IMC-REFL_DC_OUT_DQ_Omicron
ALS-X_REFL_ERR_OUT_DQ_Omicron
PEM-CS_MAG_EBAY_SUSRACK_Y_DQ_Omicron
ASC-Y_TR_B_NSUM_OUT_DQ_Omicron
ASC-AS_B_RF45_Q_PIT_OUT_DQ_Omicron
SUS-OMC_M1_ISIWIT_T_DQ_Omicron
PSL-ISS_AOM_DRIVER_MON_OUT_DQ_Omicron
LSC-PRCL_OUT_DQ_Omicron
"

for CHANNEL in $CHANNELS
do
  ls /home/detchar/triggers/ER7/$SITE/$CHANNEL/*/* | lalapps_path2cache > /home/albert.wandui/caches/$SITE/$CHANNEL.cache
done 
