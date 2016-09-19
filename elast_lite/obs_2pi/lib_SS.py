#! *** Setup rest of vars now ***
WMIN=1.400
WMAX=2.125
EXPT="e16"

Q2LIMITS={"lowQ2":[2.00,3.00],"highQ2":[3.00,5.00]}
CUMSIMS={"lowQ2":"sim4_sim5_sim6_sim7_sim8_sim13", "highQ2":"sim9_sim10_sim11_sim12"}
#! The following, SE and SECODE, should be in accordance with code in h10_2_per_non_vst_SE.py where 'SYSTEMATIC_EFFECT' are setup
#! + Therefore, should also agree with 'CUTSNCORS pretty print:' in $HOME/h10_2_per_non_vst_SE_results_logs/<Q2>_<SE>_<DATE>/h10_2_per_non_vst_SE_results.log
SES=['MM','SSBands']
SECODES={
        'MM':     [':MM-AT:gpart-pid-OFF:stat-pid-OFF:',':MM-EI:gpart-pid-OFF:stat-pid-OFF:'],
        'SSBands':[':MM-EI:gpart-pid-OFF:stat-pid-OFF:',':MM-EI:gpart-pid-ON:stat-pid-ON:']
}

#! print all constant vars
#print "*** constant vars ***"
#print "WMIN=",WMIN
#print "WMAX=",WMAX
#print "EXPT=",EXPT
#print "Q2LIMITS=",Q2LIMITS
#print "CUMSIMS=",CUMSIMS
#print "SES=",SES
#print "SECODES=",SECODES
#print "******"
