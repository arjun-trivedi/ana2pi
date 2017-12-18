#!/bin/bash

#! dirl relative to $OBSDIR_E16/SS
#dirl=(
#lowQ2_MM_091616_sim4_sim5_sim6_sim7_sim8_sim13_091916 
#lowQ2_SSBands_091616_sim4_sim5_sim6_sim7_sim8_sim13_091916
#highQ2_MM_091616_sim9_sim10_sim11_sim12_091916
#highQ2_SSBands_091616_sim9_sim10_sim11_sim12_091916
#lowQ2_cmb_non_vst_SE_091916
#highQ2_cmb_non_vst_SE_091916
#lowQ2_cmb_vst_SE_091916
#highQ2_cmb_vst_SE_091916
#)

#dirl=(
#lowQ2_MM_092516_sim4_sim5_sim6_sim7_sim8_sim13_092716 
#lowQ2_SSBands_092516_sim4_sim5_sim6_sim7_sim8_sim13_092716
#highQ2_MM_092516_sim9_sim10_sim11_sim12_092716
#highQ2_SSBands_092516_sim9_sim10_sim11_sim12_092716
#lowQ2_cmb_non_vst_SE_092716
#highQ2_cmb_non_vst_SE_092716
#lowQ2_cmb_vst_SE_092716
#highQ2_cmb_vst_SE_092716
#)

#DATE='061317'
#dirl=(
#lowQ2_MM_061217_sim4_sim5_sim6_sim7_sim8_sim13_061317
#lowQ2_SSBands_061217_sim4_sim5_sim6_sim7_sim8_sim13_061317
#highQ2_MM_061217_sim9_sim10_sim11_sim12_061317
#highQ2_SSBands_061217_sim9_sim10_sim11_sim12_061317
#lowQ2_cmb_non_vst_SE_061317
#highQ2_cmb_non_vst_SE_061317
#lowQ2_cmb_vst_SE_061317
#highQ2_cmb_vst_SE_061317
#)

#DATE='061517'
#dirl=(
#lowQ2_MM_061417_sim4_sim5_sim6_sim7_sim8_sim13_061517
#lowQ2_SSBands_061417_sim4_sim5_sim6_sim7_sim8_sim13_061517
#highQ2_MM_061417_sim9_sim10_sim11_sim12_061517
#highQ2_SSBands_061417_sim9_sim10_sim11_sim12_061517
#lowQ2_cmb_non_vst_SE_061517
#highQ2_cmb_non_vst_SE_061517
#lowQ2_cmb_vst_SE_061517
#highQ2_cmb_vst_SE_061517
#)

#DATE='071117'
#dirl=(
#lowQ2_MM_071017_sim4_sim5_sim6_sim7_sim8_sim13_071117
#lowQ2_SSBands_071017_sim4_sim5_sim6_sim7_sim8_sim13_071117
#highQ2_MM_071017_sim9_sim10_sim11_sim12_071117
#highQ2_SSBands_071017_sim9_sim10_sim11_sim12_071117
#lowQ2_cmb_non_vst_SE_071117
#highQ2_cmb_non_vst_SE_071117
#lowQ2_cmb_vst_SE_071117
#highQ2_cmb_vst_SE_071117
#)

# DATE='080317'
# dirl=(
# lowQ2_SSBands_080217_sim4_sim5_sim6_sim7_sim8_sim13_080317
# highQ2_SSBands_080217_sim9_sim10_sim11_sim12_080317
# lowQ2_cmb_non_vst_SE_080317
# highQ2_cmb_non_vst_SE_080317
# lowQ2_cmb_vst_SE_080317
# highQ2_cmb_vst_SE_080317
# )

DATE='121717'
dirl=(
lowQ2_SSBands_121517_sim4_sim5_sim6_sim7_sim8_sim13_121717
highQ2_SSBands_121517_sim9_sim10_sim11_sim12_121717
lowQ2_cmb_non_vst_SE_121717
highQ2_cmb_non_vst_SE_121717
lowQ2_cmb_vst_SE_121717
highQ2_cmb_vst_SE_121717
)

q2l=(
\[2.00,3.00\) 
\[2.00,3.00\) 
\[3.00,5.00\) 
\[3.00,5.00\)
\[2.00,3.00\)
\[3.00,5.00\)
\[2.00,3.00\)
\[3.00,5.00\)
)

printf "*** Going to plot_SE.py for following observable dirs (their corresponding q2 is also printed): ***\n"
printf "%s\n" ${dirl[@]}
printf "%s\n" ${q2l[@]}
printf "******\n"
#exit

logdir=/home/trivedia/plotSE_logs/plotSE_log_$DATE
mkdir -p $logdir
rm -rf $logdir/*

i=0
for dir in "${dirl[@]}";do
  q2=${q2l[i]}
  EXTRACT_D_E_FOR_NON_ALPHAL=("True" "False")
  for EXTRACT_D_E_FOR_NON_ALPHA in "${EXTRACT_D_E_FOR_NON_ALPHAL[@]}";do
    echo ">plot_SE.py $OBSDIR_E16/SS/$dir $q2 $EXTRACT_D_E_FOR_NON_ALPHA >& $logdir/$dir"_"$EXTRACT_D_E_FOR_NON_ALPHA.log"
    plot_SE.py $OBSDIR_E16/SS/$dir $q2 $EXTRACT_D_E_FOR_NON_ALPHA >& $logdir/$dir"_"$EXTRACT_D_E_FOR_NON_ALPHA.log
  done
  i=$((i + 1))
done

