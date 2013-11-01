#!/usr/bin/python
import sys,getopt
import re
import os,commands,subprocess,datetime,shutil,time, commands

class YieldTool():
	yielddir=h10lst=fyield=bkdir='';
	
	def __init__(self,datadir,anadir,simnum,skim,procorder,nevts,bkdirname,logfile):
		self.datadir  =datadir;
		self.anadir   =anadir;
		self.simnum   =simnum;
		self.skim     =skim;
		self.procorder=procorder;
		self.nevts    =nevts;
		self.bkdirname=bkdirname;
		self.logfile  =logfile

	def setup(self,dtype):
		if not os.path.isdir(self.anadir):
			os.makedirs(anadir)

		self.dtype=dtype;
		if self.dtype=="cooked":
			self.yielddir=os.path.join(self.anadir,self.simnum)
			self.h10lst  =os.path.join(self.yielddir,self.skim+'.lst')
			self.fyield  =os.path.join(self.yielddir,'yield.root')
			if self.bkdirname:
				self.bkdir=os.path.join(self.yielddir,self.bkdirname)
		elif self.dtype=='mc':
			self.yielddir=os.path.join(self.anadir,self.simnum)
			self.h10lst  =os.path.join(self.yielddir,self.skim+'.lst')
		 	self.fyield  =os.path.join(self.yielddir,'mcyield.root')
			if self.bkdirname:
				self.bkdir=os.path.join(self.yielddir,self.bkdirname)    	
			
					
	def display(self):
		print ('YieldTool:\nyielddir=%s\nh10lst=%s\nfyield=%s'%(self.yielddir,self.h10lst,self.fyield));
		if self.bkdir:  print ('Will update yield and backup old data in %s' % self.bkdir)
		else:           print 'No option for bkdir (i.e. if yield exists, will not overwrite)'  	 

	def makeh10lst(self):
		os.system("ls %s/%s/%s/*.root > %s"%(self.datadir,self.simnum,self.skim,self.h10lst))

	# def run(self):
	# 	if   self.dtype=='mc'    :self.run_mc()
	# 	elif self.dtype=='cooked':self.run_cooked()

	def run(self):
		print 'in make_tools::run'

		if not os.path.isdir(self.yielddir):
			os.makedirs(self.yielddir);
			self.makeh10lst();
		elif os.path.exists(self.fyield):
			if self.bkdir:
				#bkdir=os.path.join(self.yielddir,'bk__' +self.bkdir+'__'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))
				#os.makedirs(bkdir)
				#self.bkdir=os.path.join(self.yielddir,self.bkdir)
				os.makedirs(self.bkdir)
				print self.fyield,'exists! Making backup in: \n',self.bkdir
				os.system("cp -p %s/*.lst   %s" % (self.yielddir,self.bkdir))
				os.system("mv %s/*__sim.root   %s" % (self.yielddir,self.bkdir))
				os.system("mv %s/*__exp.root   %s" % (self.yielddir,self.bkdir))
				if self.dtype=='mc':
					os.system("mv %s/mcyield.root  %s" % (self.yielddir,self.bkdir))
					os.system("cp -p %s/yield.root  %s" % (self.yielddir,self.bkdir))
				elif self.dtype=='cooked':
					os.system("mv %s/yield.root  %s" % (self.yielddir,self.bkdir))
					os.system("cp -p %s/mcyield.root  %s" % (self.yielddir,self.bkdir))

				#create updated h10lst
				self.makeh10lst();
			else: 
				print self.fyield,'exists! Skipping...'
				return;
		elif not os.path.exists(self.fyield):
			if not os.path.exists(self.h10lst):
				self.makeh10lst();
		
		self.h10lst='\"'+self.h10lst+'\"'
		self.fyield='\"'+self.fyield+'\"'
		print('Going to call runSelector with following arguments')
		print('nevts=%s'%self.nevts)
		print('h10lst=%s'%self.h10lst)
		print('fyield=%s'%self.fyield)
		print('procorder=%s'%self.procorder)
		time.sleep(4)
		cmdstring='root -l -b -q \'runSelector.C(%s,%s,%s,%s); >>& %s\'' % (self.nevts,self.h10lst,self.fyield,self.procorder,self.logfile)
		print 'cmdstring = ', cmdstring;
		print '\n'
		os.system(cmdstring)