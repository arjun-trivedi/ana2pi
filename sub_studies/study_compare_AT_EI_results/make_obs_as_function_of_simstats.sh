#!/bin/bash

#! [10-16-17] 
#! + Started as a copy of $ELAST_LITE/obs_2pi/make_h10_2_per_non_vst_SE_results.sh
#! + Attempt to study how observables change with simulation statistics
#! + For SSBands-off-off (new option added to h10_2_per_non_vst_SE_results.py)

#! lowQ2
#! 1. sim4: 16.67% of sim stats
h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_SSBands_off_off_sim4                                 sim4                           SSBands-off-off &
wait
#! 2. sim4:sim5 33.33% of sim stats
h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_SSBands_off_off_sim4_sim5                            sim4:sim5                            SSBands-off-off &
wait
#! 3. sim4:sim5:sim6 50% of sim stats
h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_SSBands_off_off_sim4_sim5_sim6                       sim4:sim5:sim6                       SSBands-off-off &
wait
#! 4. sim4:sim5:sim6:sim7 66.67% of sim stats
h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_SSBands_off_off_sim4_sim5_sim6_sim7                  sim4:sim5:sim6:sim7                  SSBands-off-off &
wait
#! 5. sim4:sim5:sim6:sim7:sim8 83.33% of sim stats
h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_SSBands_off_off_sim4_sim5_sim6_sim7_sim8             sim4:sim5:sim6:sim7:sim8             SSBands-off-off &
wait
#! 6. sim4:sim5:sim6:sim7:sim8:sim13 100% of sim stats
h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_SSBands_off_off_sim4_sim5_sim6_sim7_sim8_sim13       sim4:sim5:sim6:sim7:sim8:sim13       SSBands-off-off &
wait
#! 7. sim4:sim5:sim6:sim7:sim8:sim13:sim14 100% of sim stats
h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_SSBands_off_off_sim4_sim5_sim6_sim7_sim8_sim13_sim14 sim4:sim5:sim6:sim7:sim8:sim13:sim14 SSBands-off-off &
wait

#wait

#! wait for lowQ2 to finish before moving to highQ2
wait

#! 2. highQ2
#! 1. sim9 25% of stats
h10_2_per_non_vst_SE_results.py e16 3.00 5.00 highQ2_SSBands_off_off_sim9                         sim9                         SSBands-off-off &
wait
#! sim9:sim10 50% of stats
h10_2_per_non_vst_SE_results.py e16 3.00 5.00 highQ2_SSBands_off_off_sim9_sim10                   sim9:sim10                   SSBands-off-off &
wait
#! sim9:sim10:sim11 75% of stats
h10_2_per_non_vst_SE_results.py e16 3.00 5.00 highQ2_SSBands_off_off_sim9_sim10_sim11             sim9:sim10:sim11             SSBands-off-off &
wait
#! 4. sim9:sim10:sim11:sim12 100% of stats
h10_2_per_non_vst_SE_results.py e16 3.00 5.00 highQ2_SSBands_off_off_sim9_sim10_sim11_sim12       sim9:sim10:sim11:sim12       SSBands-off-off &
wait
#! 5. sim9:sim10:sim11:sim12:sim15 100% of stats
h10_2_per_non_vst_SE_results.py e16 3.00 5.00 highQ2_SSBands_off_off_sim9_sim10_sim11_sim12_sim15 sim9:sim10:sim11:sim12:sim15 SSBands-off-off &
wait

#! wait for highQ2 to finish before exiting
#! + This will give correct estimate for time when run with the 'time' command
wait
