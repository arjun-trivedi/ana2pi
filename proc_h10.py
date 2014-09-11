#!/usr/bin/python

import sys,os,subprocess#,shlex
import getopt

"""
This file sets up all the input arguments needed by 'proc_h10' (compiled C++ program) that does the following:
h10 ----> "outdir"/"output".root
+ Note that the "outdir" should contain the h10.lst to be used by 'proc_h10'
where "output" could be
	+ d2pi
	+ deid
	+ d"substudy"

The input arguments are:
	+ h10type:<expt>:<dtyp>:<rctn>
	+ simnum
	+ output: This governs the output from 'proc_h10' and therefore sets up 'procorder','outdir' & 'fout'
	+ nentries
"""
def main(argv):
	h10type=''
	simnum=''
	output=''
	fout=''
	nentries='1000000000'
	try:
		opts, args = getopt.getopt(argv,"h",["h10type=","simnum=","output=","nentries="])
	except getopt.GetoptError:
		print('Arguments not entered correctly. Correct syntax is:')
		print('proc_h10.py --h10type=<expt>:<dtyp>:<rctn> --simnum=<simX> --output=<output> --nentries')
		sys.exit('Exiting.')
		#sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'proc_h10.py --h10type=<expt>:<dtyp>:<rctn> --simnum=<simX> --q2w=<q2wX> --use_q2w_elist=<True/False> --output=<output> --nentries'
			sys.exit()
		elif opt in ("--h10type"):
			h10type=arg
		elif opt in ("--simnum"):
			simnum=arg
		elif opt in ("--output"):
			output=arg
		elif opt in ("--nentries"):
			nentries=arg
	
	#! Make sure the input is sound
	if h10type=='' or output=='':
		sys.exit("Not all required arguments entered. Exiting.")
	expt,dtyp,rctn=h10type.split(":")
	print "expt,dtyp,rctn=",expt,dtyp,rctn
	if dtyp=='sim' and simnum=='':
		sys.exit("simnum num not entered for simulation")
	#! Print input information
	print 'h10type=%s'%h10type
	if dtyp=='sim': print 'simnum=%s'%simnum
	print 'output=%s'%output
	print 'nentries=%s'%nentries

	#! Prepare outdir,procorder
	if dtyp=='exp':
		if output=='d2pi':
			outdir=os.path.join(os.environ['D2PIDIR_EXP'])
			procorder="eid:efid:qskim:mom:pid:d2piR"
		else:
			sys.exit("output=%s not recognized"%output)
	elif dtyp=='sim':
		if output=='d2pi':
			outdir=os.path.join(os.environ['D2PIDIR_SIM'],simnum)
			procorder="d2piT:eid:efid:qskim:pid:d2piR"
		else:
			sys.exit("output=%s not recognized"%output)
	#! Prepare fout,h10lst
	h10lst=os.path.join(outdir,"h10.lst")
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	fout=os.path.join(outdir,"%s.root"%output)

	# if dtyp=='exp':
	# 	h10lst=os.path.join(os.environ['OBS_DATADIR_EXP'],"h10.lst")
	# 	outdir='$HOME/ongoing/mem_test/exp/new-h8-bng'
	# 	#outdir=os.path.join(os.environ['OBS_DATADIR'],'mem_test','exp','new-h8-bng')
	# 	# if output=='d2pi':
	# 	# 	#at-h8 outdir=os.path.join(obs_datadir,q2w)
	# 	# 	outdir=obs_datadir
	# 	# elif output=='q2welist':
	# 	# 	outdir=os.path.join(obs_datadir)
	# 	# else:
	# 	# 	sys.exit("output=%s not recognized"%output)
	# 	# h10lst=os.path.join(obs_datadir,"h10.lst")
	# 	# if use_q2w_elist:
	# 	# 	f_q2w_el=os.path.join(obs_datadir,"q2welist.root")
	# if dtyp=='sim': 
	# 	h10lst=os.path.join(os.environ['OBS_DATADIR_SIM'],'sim1',"h10.lst")
	# 	outdir='$HOME/ongoing/mem_test/sim/new-h8-bng'
	# 	#outdir=os.path.join(os.environ['OBS_DATADIR'],'mem_test','sim','new-h8-bng')
	# 	# #at-h8 obs_datadir=os.environ['OBS_DATADIR_SIM']
	# 	# obs_datadir=os.path.join(os.environ['OBS_DATADIR'],'mem_test','sim')
	# 	# if output=='d2pi':
	# 	# 	#at-h8 outdir=os.path.join(obs_datadir,simnum,q2w)
	# 	# 	outdir=obs_datadir
	# 	# elif output=='q2welist':
	# 	# 	outdir=os.path.join(obs_datadir,simnum)
	# 	# else:
	# 	# 	sys.exit("output=%s not recognized"%output)
	# 	# h10lst=os.path.join(obs_datadir,"h10.lst")
	# 	# if use_q2w_elist:
	# 	# 	f_q2w_el=os.path.join(obs_datadir,simnum,"q2welist.root")

	# if not os.path.exists(outdir):
	# 	os.makedirs(outdir)
	# fout=os.path.join(outdir,"%s.root"%output)

	#! Determine procorder

	# #! First determine proc_q2wskim
	# if q2w=='':
	# 	proc_q2wskim_dcptr=''
	# else:
	# 	num=q2w.split('q2w')[1]
	# 	proc_q2wskim_dcptr='q2wskim%s'%num	

	# if dtyp=='exp':
	# 	if output=='d2pi':
	# 		#procorder="eid:%s:efid:qskim:mom:pid:d2piR"%proc_q2wskim_dcptr
	# 		procorder="eid:efid:qskim:mom:pid:d2piR"
	# 	elif output=='q2welist':
	# 		procorder='eid:q2welist'
	# 	else:
	# 		sys.exit("output=%s not recognized"%output)
	# if dtyp=='sim':
	# 	if output=='d2pi':
	# 		#procorder="%s:d2piT:eid:efid:qskim:pid:d2piR"%proc_q2wskim_dcptr
	# 		procorder="d2piT:eid:efid:qskim:pid:d2piR"
	# 	elif output=='q2welist':
	# 		procorder='q2welist'
	# 	else:
	# 		sys.exit("output=%s not recognized"%output)


	#! Finall call proc_h10
	cmd=["proc_h10","-i",h10lst,"-t",h10type,"-p",procorder,"-o",fout,"-n",nentries]#,"-l",str(use_q2w_elist),"-m",f_q2w_el,"-q",q2w]
	print ">>>",cmd
	#subprocess.check_output(cmd)
	subprocess.call(cmd)
	
if __name__ == "__main__":
	main(sys.argv[1:])
