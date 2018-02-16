#!/usr/bin/python
import itertools,subprocess,os,sys

USAGE='run_study_yields_acceptance dbg[=False] plot_h5_stats_vst_var[=False]'

dbg="False"
if len(sys.argv)>1: #! dbg entered by user
	if    sys.argv[1]=="True":  dbg="True"
	elif  sys.argv[1]=="False": dbg="False"
	else: sys.exit('dbg=%s is not valid. usage: %s'%(sys.argv[1],USAGE))

plot_h5_stats_vst_var="False"
if len(sys.argv)>2: #!  show_rel_err_dist entered by user
	if    sys.argv[2]=="True":  plot_h5_stats_vst_var="True"
	elif  sys.argv[2]=="False": plot_h5_stats_vst_var="False"
	else: sys.exit('plot_h5_stats_vst_var=%s is not valid. usage: %s'%(sys.argv[2],USAGE))

print "dbg=",dbg
print "plot_h5_stats_vst_var=",plot_h5_stats_vst_var
#sys.exit()

show_rel_err_dist=["False","True"]

for opt_rel_err in show_rel_err_dist:
	print opt_rel_err
	exc='%s/study_yields_acceptance.py'%(os.environ['STUDY_YIELDS_ACCEPTANCE'])
	cmd=[exc,dbg,opt_rel_err,plot_h5_stats_vst_var]
	print cmd
	tool=subprocess.Popen(cmd)#,stderr=subprocess.STDOUT)
	tool.wait()

# yields_ER_commonbins=["True","False"]
# show_rel_err_dist=["True","False"]
# SSBands=['off-off','on-on']

# for item in list(itertools.product(yields_ER_commonbins,show_rel_err_dist,SSBands)):
# 	print item[0],item[1],item[2]
# 	exc='%s/study_yields_acceptance.py'%(os.environ['STUDY_YIELDS_ACCEPTANCE'])
# 	cmd=[exc,dbg,item[0],item[1],item[2]]
# 	print cmd
# 	tool=subprocess.Popen(cmd)#,stderr=subprocess.STDOUT)
# 	tool.wait()
