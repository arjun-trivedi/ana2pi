#!/usr/bin/python
import glob,os
from subprocess import Popen, PIPE, STDOUT

SUBDIR=os.path.join(os.environ['BSUBDIR'],"e1f/2pi/exp/h10_2_dlum")
if not os.path.exists(SUBDIR):
	os.makedirs(SUBDIR)
OUTDIR=os.environ['STUDY_LUM_DATADIR']

h10lstl=glob.glob(os.path.join(os.environ['HOME'],"ongoing/e1f_allruns_h10lst/*.lst"))
h10lstl.sort()
#print h10lstl
for h10lst in h10lstl:
	runnum=int(h10lst.split("/")[5].split(".")[0])
	#print runnum
	print "Going to submit job for runnum %d"%runnum
        JSUBFILE=os.path.join(SUBDIR,"%d.sub"%runnum)
	print JSUBFILE
	if os.path.exists(JSUBFILE):
		print "job.sub exists. Deleting it before creating new one."
	 	os.remove(JSUBFILE)	
	fsub=open(JSUBFILE,"w")
	fsub.write("JOBNAME: e1f_h10_2_dlum_%s\n"%runnum)
	fsub.write("PROJECT: E1F\n")
	fsub.write("SINGLE_JOB: true\n")
	fsub.write("OS: centos65\n")
	fsub.write("TRACK: analysis\n")
	fsub.write("DISK_SPACE: 500 MB\n")
	fsub.write("MEMORY: 1000 MB\n")
	fsub.write("TIME: 5\n")
	fsub.write(" \n")
        #fsub.write("COMMAND: proc_h10 -i h10.lst -t e1f:exp:2pi  -p lum -o lum.root\n"%h10lst)
        fsub.write("COMMAND: h10_2_dlum\n")
	#! Read h10lst and append $E1FD to all files
	#! and list as OTHER files for job (so that they can be got from Tape)
	fh10lst=open(h10lst)
	for fh10 in fh10lst.readlines():
		fsub.write("OTHER_FILES: %s/%s\n"%(os.environ['E1FD'],fh10.strip()))
	fsub.write("INPUT_FILES: %s\n"%h10lst)
	fsub.write("INPUT_DATA: h10.lst\n")
        fsub.write("OUTPUT_DATA: lum.root\n")
        fsub.write("OUTPUT_TEMPLATE: %s/%s.root\n"%(os.environ['STUDY_LUM_DATADIR'],runnum))

        #! Now submit the job
        cmd="jsub %s"%JSUBFILE
        #p=Popen(cmd,shell=True,stdin=PIPE,stdout=PIPE,stderr=STDOUT,close_fds=True)
	#stdout=p.stdout.read()
    	#print stdout
        print "Submitted job for runnum %d"%runnum 
