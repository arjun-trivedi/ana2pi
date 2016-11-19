#!/bin/bash

#--lQ2
$SUBSTUDIES/study_d2pi_kinematics/plot_d2pi_kin_e16.py p lQ2     >&/tmp/log_lQ2_p &
$SUBSTUDIES/study_d2pi_kinematics/plot_d2pi_kin_e16.py theta lQ2 >&/tmp/log_lQ2_theta&
$SUBSTUDIES/study_d2pi_kinematics/plot_d2pi_kin_e16.py phi lQ2   >&/tmp/log_lQ2_phi&

wait

#-- hQ2
$SUBSTUDIES/study_d2pi_kinematics/plot_d2pi_kin_e16.py p hQ2     >&/tmp/log_hQ2_p &
$SUBSTUDIES/study_d2pi_kinematics/plot_d2pi_kin_e16.py theta hQ2 >&/tmp/log_hQ2_theta &
$SUBSTUDIES/study_d2pi_kinematics/plot_d2pi_kin_e16.py phi hQ2   >&/tmp/log_hQ2_phi &
