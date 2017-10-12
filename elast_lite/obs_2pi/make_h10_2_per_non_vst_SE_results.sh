#!/bin/bash

#! 1. lowQ2
#! SSBands,MM
#h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_SSBands sim4:sim5:sim6:sim7:sim8:sim13 SSBands &
#! testdiffATEI
#h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_testdiffATEI sim4:sim5:sim6:sim7:sim8:sim13 testdiffATEI &
h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_testdiffATEI2 sim4:sim5:sim6:sim7:sim8:sim13 testdiffATEI2 &
#wait
#! [08-12-17] MM removed as a part of re-obtain-obs-3
#h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_MM      sim4:sim5:sim6:sim7:sim8:sim13 MM &
#wait

#! wait for lowQ2 to finish before moving to highQ2
wait

#! 2. highQ2
#! SSBands,MM
#h10_2_per_non_vst_SE_results.py e16 3.00 5.00 highQ2_SSBands sim9:sim10:sim11:sim12 SSBands &
#! testdiffATEI
#h10_2_per_non_vst_SE_results.py e16 3.00 5.00 highQ2_testdiffATEI sim9:sim10:sim11:sim12 testdiffATEI &
#h10_2_per_non_vst_SE_results.py e16 3.00 5.00 highQ2_testdiffATEI2 sim9:sim10:sim11:sim12 testdiffATEI2 &
#wait
#! [08-12-17] MM removed as a part of re-obtain-obs-3
#h10_2_per_non_vst_SE_results.py e16 3.00 5.00 highQ2_MM      sim9:sim10:sim11:sim12 MM &
#wait

#! wait for highQ2 to finish before exiting
#! + This will give correct estimate for time when run with the 'time' command
wait
