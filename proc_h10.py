#!/usr/bin/python

import sys,os,getopt

def main(argv):
	h10type=''
	q2w=''
	output=''
	try:
		opts, args = getopt.getopt(argv,"h",["h10type=","q2w=","output="])
	except getopt.GetoptError:
		print 'proc_h10.py h10type= q2w= output='
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'proc_h10.py --h10type=<expt>:<dtyp>:<rctn>, --q2w=<q2wX>, --output=<output>'
			sys.exit()
		elif opt in ("--h10type"):
			h10type=arg
		elif opt in ("--q2w"):
			q2w=arg
		elif opt in ("--output"):
			output=arg
	if h10type=='' or q2w=='' or output=='':
		sys.exit("Not all required arguments entered. Exiting.")
	print 'h10type=%s'%h10type
	print 'q2w=%s'%q2w
	print 'output=%s'%output

	expt,dtyp,rctn=h10type.split(":")
	print "expt,dtyp,rctn=",expt,dtyp,rctn

	host=os.environ['HOST']
	if "gothe14" in host:
		print "HOST=%s"%host
	elif "ifarm" in host:
		print "HOST=%s"%host
	else:
		sys.exit("HOST not determined. Exiting")
	#exec="proc_h10 -i %s/h10.lst"

if __name__ == "__main__":
	main(sys.argv[1:])