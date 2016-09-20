#!/bin/bash

#! dirl relative to $OBSDIR_E16/SS
dirl=(
lowQ2_MM_091616_sim4_sim5_sim6_sim7_sim8_sim13_091916 
lowQ2_SSBands_091616_sim4_sim5_sim6_sim7_sim8_sim13_091916
highQ2_MM_091616_sim9_sim10_sim11_sim12_091916
highQ2_SSBands_091616_sim9_sim10_sim11_sim12_091916
lowQ2_cmb_non_vst_SE_091916
highQ2_cmb_non_vst_SE_091916
lowQ2_cmb_vst_SE_091916
highQ2_cmb_vst_SE_091916
)
printf "*** Going to plot_SE.py for following observable dirs: ***"
printf "%s\n" ${dirl[@]}
printf "******\n"

logdir=/home/trivedia/plot_SE_logs
mkdir -p $logdir
rm -rf $logdir/*
for dir in "${dirl[@]}";do
  EXTRACT_D_E_FOR_NON_ALPHAL=("True" "False")
  for EXTRACT_D_E_FOR_NON_ALPHA in "${EXTRACT_D_E_FOR_NON_ALPHAL[@]}";do
    echo ">plot_SE.py $OBSDIR_E16/SS/$dir $EXTRACT_D_E_FOR_NON_ALPHA >& $logdir/$dir"_"$EXTRACT_D_E_FOR_NON_ALPHA.log"
    plot_SE.py $OBSDIR_E16/SS/$dir $EXTRACT_D_E_FOR_NON_ALPHA >& $logdir/$dir"_"$EXTRACT_D_E_FOR_NON_ALPHA.log
  done
done

