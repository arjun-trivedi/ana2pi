#!/usr/bin/python

import sys,os,subprocess#,shlex
import getopt

def main(argv):
	h10type=''
	sim_num=''
	q2w=''
	use_q2w_elist=False
	f_q2w_el=''
	output=''
	fout=''
	nentries='1000000000'
	try:
		opts, args = getopt.getopt(argv,"h",["h10type=","sim_num=","q2w=","output=","use_q2w_elist=","nentries="])
	except getopt.GetoptError:
		print('Arguments not entered correctly. Correct syntax is:')
		print('proc_h10.py --h10type=<expt>:<dtyp>:<rctn> --sim_num=<simX> --q2w=<q2wX> --use_q2w_elist=<True/False> --output=<output> --nentries')
		sys.exit('Exiting.')
		#sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'proc_h10.py --h10type=<expt>:<dtyp>:<rctn> --sim_num=<simX> --q2w=<q2wX> --use_q2w_elist=<True/False> --output=<output> --nentries'
			sys.exit()
		elif opt in ("--h10type"):
			h10type=arg
		elif opt in ("--sim_num"):
			sim_num=arg
		elif opt in ("--q2w"):
			q2w=arg
		elif opt in ("--use_q2w_elist"):
			use_q2w_elist=arg
		elif opt in ("--output"):
			output=arg
		elif opt in ("--nentries"):
			nentries=arg
	
	#! Print all the input data
	print 'h10type=%s'%h10type
	print 'sim_num=%s'%sim_num
	print 'q2w=%s'%q2w
	print 'use_q2w_elist=%s'%use_q2w_elist
	print 'output=%s'%output
	print 'nentries=%s'%nentries

	if h10type=='' or output=='':
		sys.exit("Not all required arguments entered. Exiting.")

	expt,dtyp,rctn=h10type.split(":")
	print "expt,dtyp,rctn=",expt,dtyp,rctn
	if dtyp=='sim' and sim_num=='':
		sys.exit("sim_num num not entered for simulation")

	#! Prepare h10lst,output,f_q2w_el
	if dtyp=='exp':
		obs_datadir=os.environ['OBS_DATADIR_EXP']
		if output=='d2pi':
			outdir=os.path.join(obs_datadir,q2w)
		elif output=='q2welist':
			outdir=os.path.join(obs_datadir)
		else:
			sys.exit("output=%s not recognized"%output)
		h10lst=os.path.join(obs_datadir,"h10.lst")
		if use_q2w_elist:
			f_q2w_el=os.path.join(obs_datadir,"q2welist.root")
	if dtyp=='sim': 
		obs_datadir=os.environ['OBS_DATADIR_SIM']
		if output=='d2pi':
			outdir=os.path.join(obs_datadir,sim_num,q2w)
		elif output=='q2welist':
			outdir=os.path.join(obs_datadir,sim_num)
		else:
			sys.exit("output=%s not recognized"%output)
		h10lst=os.path.join(obs_datadir,sim_num,"h10.lst")
		if use_q2w_elist:
			f_q2w_el=os.path.join(obs_datadir,sim_num,"q2welist.root")

	if not os.path.exists(outdir):
		os.makedirs(outdir)
	fout=os.path.join(outdir,"%s.root"%output)

	#! Determine procorder

	#! First determine proc_q2wskim
	if q2w=='':
		proc_q2wskim_dcptr=''
	else:
		num=q2w.split('q2w')[1]
		proc_q2wskim_dcptr='q2wskim%s'%num	

	if dtyp=='exp':
		if output=='d2pi':
			#procorder="eid:%s:efid:qskim:mom:pid:d2piR"%proc_q2wskim_dcptr
			procorder="eid:efid:qskim:mom:pid:d2piR"
		elif output=='q2welist':
			procorder='eid:q2welist'
		else:
			sys.exit("output=%s not recognized"%output)
	if dtyp=='sim':
		if output=='d2pi':
			procorder="%s:d2piT:eid:efid:qskim:pid:d2piR"%proc_q2wskim_dcptr
			#procorder="d2piT:eid:efid:qskim:pid:d2piR"
		elif output=='q2welist':
			procorder='q2welist'
		else:
			sys.exit("output=%s not recognized"%output)


	#! Finall call proc_h10
	cmd=["proc_h10","-i",h10lst,"-t",h10type,"-p",procorder,"-o",fout,"-l",str(use_q2w_elist),"-m",f_q2w_el,"-q",q2w,"-n",nentries]
	print ">>>",cmd
	subprocess.check_output(cmd)
	
if __name__ == "__main__":
	main(sys.argv[1:])