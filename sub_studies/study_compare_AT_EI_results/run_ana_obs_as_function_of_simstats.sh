#!/bin/bash

#! run with SIMSTATS_SHOW_REL_ERR_DIST==False
$SUBSTUDIES/study_compare_AT_EI_results/ana_obs_as_function_of_simstats.py False False False
wait
#! run with SIMSTATS_SHOW_REL_ERR_DIST==True
$SUBSTUDIES/study_compare_AT_EI_results/ana_obs_as_function_of_simstats.py False True False
