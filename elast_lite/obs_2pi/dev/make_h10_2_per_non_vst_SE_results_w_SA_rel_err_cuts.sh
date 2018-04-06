#!/bin/bash

#! Define SA-rel-err-cut
CUT_VALS=(0.50 0.65 0.95)
CUT_NAMES=(050  065  095)

#! Set up some constants
SKIP_MAKE_D2PI=True
SKIP_DOBS_R2=False
ACC_CUT=-1

icut=0
for cutval in "${CUT_VALS[@]}";do
	cutname=${CUT_NAMES[icut]}
	echo "cutval,cutname,SKIP_MAKE_D2PI,SKIP_DOBS_R2,ACC_CUT="$cutval,$cutname,$SKIP_MAKE_D2PI,$SKIP_DOBS_R2,$ACC_CUT
	#echo lowQ2_SSBands_SA_rel_err_cut_"$cutname"_032718

	#! 1. lowQ2
	#! SSBands,MM
	h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_SSBands_SA_rel_err_cut_"$cutname" sim4:sim5:sim6:sim7:sim8:sim13:sim14 SSBands 032718 $SKIP_MAKE_D2PI $SKIP_DOBS_R2 $ACC_CUT $cutval&
	#wait
	#! [08-12-17] MM removed as a part of re-obtain-obs-3
	#h10_2_per_non_vst_SE_results.py e16 2.00 3.00 lowQ2_MM      sim4:sim5:sim6:sim7:sim8:sim13 MM &
	#wait

	#! wait for lowQ2 to finish before moving to highQ2
	wait

	#! 2. highQ2
	#! SSBands,MM
	h10_2_per_non_vst_SE_results.py e16 3.00 5.00 highQ2_SSBands_SA_rel_err_cut_"$cutname" sim9:sim10:sim11:sim12:sim15 SSBands 032718 $SKIP_MAKE_D2PI $SKIP_DOBS_R2 $ACC_CUT $cutval&
	#wait
	#! [08-12-17] MM removed as a part of re-obtain-obs-3
	#h10_2_per_non_vst_SE_results.py e16 3.00 5.00 highQ2_MM      sim9:sim10:sim11:sim12 MM &
	#wait

	#! wait for highQ2 to finish before exiting
	#! + This will give correct estimate for time when run with the 'time' command
	wait

	icut=$((icut + 1))
done
