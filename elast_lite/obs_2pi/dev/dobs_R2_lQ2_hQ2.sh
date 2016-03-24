#!/bin/bash
DBG=${1:-False}

viewl=(EC_SF_ST EC_EF)

for view in "${viewl[@]}";do
  dobs_R2 mthd1:mthd2:True $OBSDIR_E16/Archive/obs_dflt_eff_scpd_at_mod_dflt_lowQ2_SS3_022216/ sim4_sim5_sim6_sim7 $view 2.00 3.00 1.400 2.125 e16 $DBG >& /tmp/lQ2_$view.log

  dobs_R2 mthd1:mthd2:True $OBSDIR_E16/Archive/obs_dflt_eff_scpd_at_mod_dflt_highQ2_SS3_022216/ sim9 $view 3.00 5.00 1.400 2.125 e16 $DBG >& /tmp/hQ2_$view.log
done
