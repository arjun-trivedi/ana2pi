#!/bin/bash
#! $1=DATE when h10_2_per_non_vst_SE_results were obtain
#! $2=DBG. False by default
DATE=$1
DBG=${2:-"False"}

if [[ -z "$DATE" ]]; then 
  echo "No date specified"
  exit
fi

echo "DATE="$DATE
echo "DBG="$DBG

#! Setup 'TODAY_DATE' which will be used to tag SS results from 1.,2., and 3.
TODAY_DATE=`date +%m%d%y`

#! Prepare logdir
if [[ "$DBG" == "True" ]]; then
  LOGDIR="/home/trivedia/SSlogs/dbg/SSlog_$TODAY_DATE"
elif [[ "$DBG" == "False" ]]; then
  LOGDIR="/home/trivedia/SSlogs/SSlog_$TODAY_DATE"
else
  echo "DBG="$DBG is not valid
  exit
fi
mkdir -p $LOGDIR
rm -rf $LOGDIR/*
#exit

#! 1. Obtain per_non_vst_SE results
logdir=$LOGDIR/dobs_per_non_vst_SE_result
mkdir -p $logdir
rm -rf $logdir/*
dobs_per_non_vst_SE_result.py $DATE lowQ2  MM      $DBG >& $logdir/lowQ2_MM.log 
dobs_per_non_vst_SE_result.py $DATE lowQ2  SSBands $DBG >& $logdir/lowQ2_SSBands.log 
dobs_per_non_vst_SE_result.py $DATE highQ2 MM      $DBG >& $logdir/highQ2_MM.log 
dobs_per_non_vst_SE_result.py $DATE highQ2 SSBands $DBG >& $logdir/highQ2_SSBands.log 

#! Wait for processes in 1 to finish before moving on
#wait

#! 2. Using per_non_vst_SE results from 1, obtain cmb_non_vst_SE
logdir=$LOGDIR/dobs_cmb_non_vst_SE_result
mkdir -p $logdir
rm -rf $logdir/*
dobs_cmb_non_vst_SE_result.py $DATE $TODAY_DATE lowQ2  $DBG >& $logdir/lowQ2.log
dobs_cmb_non_vst_SE_result.py $DATE $TODAY_DATE highQ2 $DBG >& $logdir/highQ2.log

#! Wait for processes in 2 to finish before moving on
#wait


#! 3. Using cmb_non_vst_SE from 1 to obtain cmb_vst_SE
logdir=$LOGDIR/dobs_cmb_vst_SE_result
mkdir -p $logdir
rm -rf $logdir/*
dobs_cmb_vst_SE_result.py $TODAY_DATE lowQ2  $DBG >& $logdir/lowQ2.log
dobs_cmb_vst_SE_result.py $TODAY_DATE highQ2 $DBG >& $logdir/highQ2.log

#! Wait for processes in 3 to finish before moving on
#wait

#! To make SE plots for a particular OBSDIR
#! NOTE that the date has to be manually adjusted
#EXTRACT_D_E_FOR_NON_ALPHAL=("True" "False")
#for EXTRACT_D_E_FOR_NON_ALPHA in "${EXTRACT_D_E_FOR_NON_ALPHAL[@]}";do
  #! 1. First do each of the per_non_vst_SE
#  plot_SE.py $OBSDIR_E16/SS/dbg/lowQ2_fixMM2cut_062416_sim4_sim5_sim6_sim7_sim8_sim13_$TODAY_DATE/ $EXTRACT_D_E_FOR_NON_ALPHA >& $logdir/plot_SE_per_non_vst_SE_result_SSBands.log &
#  plot_SE.py $OBSDIR_E16/SS/dbg/lowQ2_CC_cut_080216_sim4_sim5_sim6_sim7_sim8_sim13_$TODAY_DATE/ $EXTRACT_D_E_FOR_NON_ALPHA >& $logdir/plot_SE_per_non_vst_SE_result_CCcut.log &
  #! 2. Now do the cmb_non_vst_SE
#  plot_SE.py $OBSDIR_E16/SS/dbg/lowQ2_cmb_non_vst_SE_$TODAY_DATE $EXTRACT_D_E_FOR_NON_ALPHA >& $logdir/plot_SE_cmb_non_vst_SE_result.log&
  #! 3. Now do the cmb_vst_SE
#  plot_SE.py $OBSDIR_E16/SS/dbg/lowQ2_cmb_vst_SE_$TODAY_DATE $EXTRACT_D_E_FOR_NON_ALPHA >& $logdir/plot_SE_cmb_vst_SE_result.log&
#done
