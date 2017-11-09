#!/bin/bash

#! Following are dummy values and not thought through yet, but are needed as required inputs
expctd_diff=1
expctd_diff_tlrnc=1

#! Remove old logs for {lowQ2,highQ2}*{0.005,0.01}
rm /tmp/lowQ2_acc005.log
rm /tmp/lowQ2_acc01.log
rm /tmp/highQ2_acc005.log
rm /tmp/highQ2_acc01.log

#! lowQ2*{0.005,0.01}
./plot_diff_obs_as_function_of_acc_cut.py $OBSDIR_E16/lowQ2_SSBands_off_off_sim4_sim5_sim6_sim7_sim8_sim13_102617/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13 $OBSDIR_E16/lowQ2_SSBands_off_off_sim4_sim5_sim6_sim7_sim8_sim13_102617_acc_cut_005/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13 $expctd_diff $expctd_diff_tlrnc >& /tmp/lowQ2_acc005.log&
./plot_diff_obs_as_function_of_acc_cut.py $OBSDIR_E16/lowQ2_SSBands_off_off_sim4_sim5_sim6_sim7_sim8_sim13_102617/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13 $OBSDIR_E16/lowQ2_SSBands_off_off_sim4_sim5_sim6_sim7_sim8_sim13_102617_acc_cut_01/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13 $expctd_diff $expctd_diff_tlrnc  >& /tmp/lowQ2_acc01.log&

#! highQ2*{0.005,0.01}
./plot_diff_obs_as_function_of_acc_cut.py $OBSDIR_E16/highQ2_SSBands_off_off_sim9_sim10_sim11_sim12_102717/cutsncors1/sim9_sim10_sim11_sim12 $OBSDIR_E16/highQ2_SSBands_off_off_sim9_sim10_sim11_sim12_102717_acc_cut_005/cutsncors1/sim9_sim10_sim11_sim12 $expctd_diff $expctd_diff_tlrnc >& /tmp/highQ2_acc005.log&
./plot_diff_obs_as_function_of_acc_cut.py $OBSDIR_E16/highQ2_SSBands_off_off_sim9_sim10_sim11_sim12_102717/cutsncors1/sim9_sim10_sim11_sim12 $OBSDIR_E16/highQ2_SSBands_off_off_sim9_sim10_sim11_sim12_102717_acc_cut_01/cutsncors1/sim9_sim10_sim11_sim12 $expctd_diff $expctd_diff_tlrnc >& /tmp/highQ2_acc01.log&
