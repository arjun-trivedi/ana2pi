#!/usr/bin/python

import sys,os,subprocess#,shlex
import getopt

"""
The job of this program is to set up the input arguments and call the underlying compiled C++ program
'proc_h10' that processes a chain of h10s to produce an output as illustrated below:
h10 -> proc1 -> proc2 -> ... -> procX -> fout

'proc_h10' is called as:
>proc_h10 -i h10lst -t h10type -p procorder -o fout -n nentries
where
	+ h10lst = TChain of h10s
	+ h10type = Tells 'proc_h10' about the type of h10s in the chain;<expt>:<dtyp>:<rctn>. Among other possible
	uses, this information is chiefly used to appropriately Bind the TBranches, which vary with h10type,
	to C++ Arrays created for them by 'proc_h10'
	+ procorder = Sequence of EProcessors that process the data in the h10 chain
	+ fout = name of output file created by 'proc_h10'
	+ nentries = total number of entries in the chain to process

The main reason for having this program is to efficiently call 'proc_h10' based on the different ways
in which I would like to process h10s. The "different ways in which I would like to process h10s" is encoded
in the desired <output> from proc_h10, which could be:
	+ d2piR,d2piT
	+ deid
	+ d<substudy>
Depending on the desired <output>, this program sets up the necessary inputs to 'proc_h10': 
'h10lst', procorder'& 'fout'(=<outdir>/<output>.root). The sytax to call this program is:
proc_h10.py --h10type=<expt>:<dtyp>:<rctn> --simnum=<simX> --output=<output> --nentries --debug=[false]

+ The input arguments are:
	+ h10type = <expt>:<dtyp>:<rctn>; directly passed to proc_h10
	+ output = This is the main input which is used to set up the call to 'proc_h10':
	'h10lst','procorder' & 'fout'(=<outdir>/<output>.root)
	+ simnum = This is mainly used to set up the appropriate <outdir> depending on the simX
	+ h10lst=optional. If not provided, as is the case for d2piR,d2piT, h10lst is expected to be in <outdir>
	+ debug = This creates output in <outdir>/debug
	+ nentries = optional argument that can be used for debugging

+ Note that for d2piR,d2piT, <outdir> should contain the h10.lst to be used by 'proc_h10'

+ Currently the following outputs are being made from h10:
	+ d2piR: for exp (h10_2_d2piR-exp)
	+ d2piR,d2piT: for sim (h10_2_d2piR-sim,h10_2_d2piT-sim)
	+ h10skim: for exp (h10_2_h10skim-exp)
	+ memtest
"""
def main(argv):
	h10type=''
	simnum=''
	output=''
	h10lst=''
	debug='false'
	fout=''
	nentries='1000000000'
	
	try:
		opts, args = getopt.getopt(argv,"h",["h10type=","simnum=","output=", "h10lst=", "nentries=","debug="])
	except getopt.GetoptError:
		print('Arguments not entered correctly. Correct syntax is:')
		print('proc_h10.py --h10type=<expt>:<dtyp>:<rctn> --simnum=<simX> --output=<output> --h10lst=<h10lst> --nentries --debug=[false]')
		sys.exit('Exiting.')
		#sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'proc_h10.py --h10type=<expt>:<dtyp>:<rctn> --simnum=<simX> --output=<output> --h10lst=<h10lst> --nentries --debug=[false]'
			sys.exit()
		elif opt in ("--h10type"):
			h10type=arg
		elif opt in ("--simnum"):
			simnum=arg
		elif opt in ("--output"):
			output=arg
		elif opt in ("--h10lst"):
			h10lst=arg
		elif opt in ("--nentries"):
			nentries=arg
		elif opt in ("--debug"):
			debug=arg
	
	#! Make sure the input arguments are sound
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
	print 'debug=%s'%debug
	if output=='h10skim': 
		if h10lst=='':
			print "For output=h10skim, input argument h10lst is needed!"
			sys.exit()
		else: print 'h10lst=%s'%h10lst

	#! Prepare h10lst,procorder,fout
	if dtyp=='exp':
		if output=='d2piR':
			if debug=='true':outdir=os.path.join(os.environ['D2PIDIR_EXP'],'debug')
			else:            outdir=os.path.join(os.environ['D2PIDIR_EXP'])	
			if not os.path.exists(outdir):
				print "%s does not exist! Please create it and put appropriate h10.lst in it"%outdir
				sys.exit()
			h10lst=os.path.join(outdir,"h10.lst")
			procorder="eid:efid:qskim:mom:pid:d2piR"
			fout=os.path.join(outdir,"%s.root"%output)
		elif output=='h10skim':
			# if debug=='true':outdir=os.path.join(os.environ['H10SKIMDIR_EXP'],debug)
			# else:            outdir=os.path.join(os.environ['H10SKIMDIR_EXP'])
			#! Since these jobs will be run on the farm, the jobs local dir will store the output
			#! and the submit script will appropriately rename the output and put it in H10SKIMDIR_EXP
			outdir=("./") 
			#! h10lst is input by user
			procorder='skimh10'
			fout=os.path.join(outdir,"%s.root"%output)
		elif output=='memtest':
			outdir=os.path.join(os.environ['D2PIDIR'],'memtest')
			h10lst=os.path.join(outdir,"h10.lst")
			procorder="eid:efid:qskim:mom:pid:d2piR"
			fout=os.path.join(outdir,"%s.root"%output)
		else:
			sys.exit("output=%s not recognized"%output)
	elif dtyp=='sim':
		if output.find('d2pi')>=0:
			if debug=='true':outdir=os.path.join(os.environ['D2PIDIR_SIM'],simnum,'debug')
			else:            outdir=os.path.join(os.environ['D2PIDIR_SIM'],simnum)
			if not os.path.exists(outdir):
				print "%s does not exist! Please create it and put appropriate h10.lst in it"%outdir
				sys.exit()
			h10lst=os.path.join(outdir,"h10.lst")
			if   output.find("T")>=0:  procorder="d2piT" 
			elif output.find("R")>=0:  procorder="eid:efid:qskim:pid:d2piR" 
			else:                      sys.exit("output=%s is invalid. output can only be d2piR or d2piT!"%output)
			fout=os.path.join(outdir,"%s.root"%output)
		elif output=='memtest':
			outdir=os.path.join(os.environ['D2PIDIR'],'memtest')
			h10lst=os.path.join(outdir,"h10.lst")
			procorder="eid:efid:qskim:pid:d2piR"
			fout=os.path.join(outdir,"%s.root"%output)
		else:
			sys.exit("output=%s not recognized"%output)

	#! Prepare fout,h10lst
	# if not os.path.exists(outdir):
	# 	print "%s does not exist! Please create it"%outdir
	# 	if output.find('d2pi')>=0:
	# 		print "And also put appropriate h10.lst in it!"
	# 	sys.exit()
	# if output.find('d2pi')>=0:
	# 	h10lst=os.path.join(outdir,"h10.lst")
	# fout=os.path.join(outdir,"%s.root"%output)

	#! Finall call proc_h10
	cmd=["proc_h10","-i",h10lst,"-t",h10type,"-p",procorder,"-o",fout,"-n",nentries]#,"-l",str(use_q2w_elist),"-m",f_q2w_el,"-q",q2w]
	print ">>>",cmd
	#subprocess.check_output(cmd)
	subprocess.call(cmd)
	
if __name__ == "__main__":
	main(sys.argv[1:])
