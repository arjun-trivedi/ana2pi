#!/bin/bash

#! Define SA-rel-err-cut
CUT_VALS=(0.50 0.65 0.95)
CUT_NAMES=(050  065  095)

#! cut observables: $OBSDIR_E16/lowQ2(highQ2)_SSBands_SA_rel_err_cut_"$cutname"_032718

#! 1. Plot diff of 're-obtain-obs-5' (current official observables: lowQ2(highQ2)_SSBands_122717) with obs-with-cuts
#! + Note that a direct comparison for with and witout SA-rel-err cut is valid only with respect to 're-obtain-obs-6' because these cut-observables obtained by applying these cuts on 're-obtain-obs-6' data
re_obtain_obs_5_date=122717
for cutname in "${CUT_NAMES[@]}";do
	$SUBSTUDIES/study_obs_as_function_of_SA_rel_err_cut/plot_diff_obs_SA_rel_err_cuts.py $OBSDIR_E16/lowQ2_SSBands_"$re_obtain_obs_5_date"/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13 $OBSDIR_E16/highQ2_SSBands_"$re_obtain_obs_5_date"/cutsncors1/sim9_sim10_sim11_sim12 $OBSDIR_E16/lowQ2_SSBands_SA_rel_err_cut_"$cutname"_032718/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13_sim14 $OBSDIR_E16/highQ2_SSBands_SA_rel_err_cut_"$cutname"_032718/cutsncors1/sim9_sim10_sim11_sim12_sim15 0 0 re_obtain_obs_5_diff_cuts__SSBands_122717_diff_SSBands_SA_rel_err_cut_"$cutname"_032718
done

#! Plot diff of 're-obtain-obs-6' (added sim14(lQ2), sim15(hQ2): lowQ2(highQ2)_SSBands_032718) with obs-with-cuts
re_obtain_obs_6_date=032718
for cutname in "${CUT_NAMES[@]}";do
	$SUBSTUDIES/study_obs_as_function_of_SA_rel_err_cut/plot_diff_obs_SA_rel_err_cuts.py $OBSDIR_E16/lowQ2_SSBands_"$re_obtain_obs_6_date"/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13_sim14 $OBSDIR_E16/highQ2_SSBands_"$re_obtain_obs_6_date"/cutsncors1/sim9_sim10_sim11_sim12_sim15 $OBSDIR_E16/lowQ2_SSBands_SA_rel_err_cut_"$cutname"_032718/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13_sim14 $OBSDIR_E16/highQ2_SSBands_SA_rel_err_cut_"$cutname"_032718/cutsncors1/sim9_sim10_sim11_sim12_sim15 0 0 re_obtain_obs_6_diff_cuts__SSBands_032718_diff_SSBands_SA_rel_err_cut_"$cutname"_032718
done

#! 3. Plot diff between various cut levels
#! Use preferred cut of 0.65 and compare 0.95 and 0.50 to it
for cutname in "${CUT_NAMES[@]}";do
	if [ "$cutname" -eq 065 ]; then
		#echo "$cutname continuing"
		continue
	fi
	#echo "$cutname processing"
	$SUBSTUDIES/study_obs_as_function_of_SA_rel_err_cut/plot_diff_obs_SA_rel_err_cuts.py $OBSDIR_E16/lowQ2_SSBands_SA_rel_err_cut_065_032718/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13_sim14 $OBSDIR_E16/highQ2_SSBands_SA_rel_err_cut_065_032718/cutsncors1/sim9_sim10_sim11_sim12_sim15 $OBSDIR_E16/lowQ2_SSBands_SA_rel_err_cut_"$cutname"_032718/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13_sim14 $OBSDIR_E16/highQ2_SSBands_SA_rel_err_cut_"$cutname"_032718/cutsncors1/sim9_sim10_sim11_sim12_sim15 0 0 SSBands_SA_rel_err_cut_065_032718_diff_SSBands_SA_rel_err_cut_"$cutname"_032718
done
