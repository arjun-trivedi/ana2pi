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
			#procorder="d2piT:eid:efid:qskim:pid:d2piR"
			procorder="d2piT"
		else:
			sys.exit("output=%s not recognized"%output)
	#! Prepare fout,h10lst
	h10lst=os.path.join(outdir,"h10.lst")
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	fout=os.path.join(outdir,"%s.root"%output)

	#! Finall call proc_h10
	cmd=["proc_h10","-i",h10lst,"-t",h10type,"-p",procorder,"-o",fout,"-n",nentries]#,"-l",str(use_q2w_elist),"-m",f_q2w_el,"-q",q2w]
	print ">>>",cmd
	#subprocess.check_output(cmd)
	subprocess.call(cmd)
	
if __name__ == "__main__":
	main(sys.argv[1:])
