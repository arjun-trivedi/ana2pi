#!/usr/bin/python
import itertools,subprocess,os,sys

USAGE='run_study_yields_acceptance dbg[=False]'

dbg="False"
if len(sys.argv)>1: #! dbg entered by user
	if    sys.argv[1]=="True":  dbg="True"
	elif  sys.argv[1]=="False": dbg="False"
	else: sys.exit('dbg=%s is not valid. usage: %s'%(sys.argv[1],USAGE))
print "dbg=",dbg
#sys.exit()

yields_ER_commonbins=["True","False"]
show_rel_err_dist=["True","False"]
SSBands=['off-off','on-on']

for item in list(itertools.product(yields_ER_commonbins,show_rel_err_dist,SSBands)):
	print item[0],item[1],item[2]
	exc='%s/study_yields_acceptance.py'%(os.environ['STUDY_YIELDS_ACCEPTANCE'])
	cmd=[exc,dbg,item[0],item[1],item[2]]
	print cmd
	tool=subprocess.Popen(cmd)#,stderr=subprocess.STDOUT)
	tool.wait()
