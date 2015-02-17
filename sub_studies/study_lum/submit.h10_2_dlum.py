#!/usr/bin/python
import glob,os
from subprocess import Popen, PIPE, STDOUT

SUBDIR="/tmp/e1f/2pi/exp/h10_2_dlum"#${BSUBDIR}/e1f/2pi/exp/h10_2_dlum
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
	ff=open(JSUBFILE,"w")
	ff.write("JOBNAME: e1f_h10_2_dlum_%s\n"%runnum)
	ff.write("PROJECT: E1F\n")
	ff.write("SINGLE_JOB: true\n")
	ff.write("OS: centos65\n")
	ff.write("TRACK: analysis\n")
	ff.write("DISK_SPACE: 500 MB\n")
	ff.write("MEMORY: 1000 MB\n")
	ff.write("TIME: 2\n")
	ff.write(" \n")
	#echo "COMMAND: h10_2_h10skim-exp "$h10lst  >> ${JSUBFILE}
	#echo "COMMAND: h10_2_h10skim-exp /tmp/atrivedi/$idtfr.lst"  >> ${JSUBFILE} 
	ff.write("COMMAND: proc_h10 -i %s -t e1f:exp:2pi  -p lum -o %s/%d.root\n"%(h10lst,os.environ['STUDY_LUM_DATADIR'],runnum))
        ff.write("OTHER_FILES: %s\n"%(os.path.join(os.environ['ANA2PI'],"proc_h10")))
	#read h10lst and list as OTHER files from E1FD
	fh10lst=open(h10lst)
	for rf in fh10lst.readlines():
		ff.write("OTHER_FILES: %s/%s\n"%(os.environ['E1FD'],rf.strip()))
	ff.write("INPUT_FILES: %s\n"%h10lst)
	ff.write("INPUT_DATA: h10.lst\n")
