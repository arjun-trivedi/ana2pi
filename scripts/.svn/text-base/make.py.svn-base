#!/usr/bin/python
import sys, getopt
import re
import os, commands, subprocess, datetime, shutil,time

import make_tools

class MyOutput():
   def __init__(self, logfile):
      self.stdout = sys.stdout
      self.log = open(logfile, 'w')

   def write(self, text):
      self.stdout.write(text)
      self.log.write(text)
      self.log.flush()
 
   def close(self):
      self.stdout.close()
      self.log.close()

#some global variables

make_yield=make_mcyield=make_h10skim=make_h10skim_2piCuts=make_q2wsel=False;
h10type=''; #e1fd1/e1fd2/e1fs1/e1fs2
skim='h10';
sim=exp=False;
bkdirname=''
#following passed to runSelector.C
nevts= 1000000000;
datadir=anadir=procorder=''

def parseopts(argv):
   global make_yield,make_mcyield,make_h10skim,make_h10skim_2piCuts,make_q2wsel
   global h10type;
   global skim;
   global sim,exp;
   global bkdirname;
   global nevts;
   global datadir,anadir,procorder;
   
   try:
      opts, args = getopt.getopt(sys.argv[1:],"h",["yield=", "mcyield=", "h10skim=", "h10skim_2piCuts", 
                                                   "nevts=", "skim=", "bkdirname="])
   except getopt.GetoptError:
      print 'Invalid arguments given to make';
      print 'make.py --yield/mcyield/h10skim/h10skim_2piCuts/q2wsel --skim --nevts --bkdirname'
      sys.exit(2)
   for opt, arg in opts:
      if opt=='-h':
         print 'make.py --yield/mcyield/h10skim/h10skim_2piCuts/q2wsel --skim --nevts --bkdirname'
         sys.exit()
      elif opt in ("--yield"):
         make_yield=True;
         h10type=arg; 
      elif opt in ("--mcyield"):
         make_mcyield=True;
         h10type=arg;
      elif opt in ("--h10skim"):
         make_h10skim=True;
         h10type=arg;
      elif opt in ("--h10skim_2piCuts"):
         make_h10skim_2piCuts=True;
         h10type=arg;
      elif opt in ("--q2wsel"):
         make_q2wsel=True;
         h10type=arg;
      elif opt in ("--skim"):
         skim=arg;
      elif opt in ("--nevts"):
         nevts=arg;
      elif opt in ("--bkdirname"):
         bkdirname=arg
         #bkdirname=(arg+'__'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))   

def init():
   global make_yield,make_mcyield,make_h10skim,make_h10skim_2piCuts,make_q2wsel
   global h10type;
   global skim;
   global sim,exp;
   global bkdirname;
   global nevts;
   global datadir,anadir,procorder;

   #check to see if h10type is valid
   if (h10type!="e1fd1")&(h10type!="e1fd2")&(h10type!="e1fs1")&(h10type!="e1fs2"):
      print('h10type %s not recognized')%h10type;
      print'Valid h10type = e1fd1/e1fd2/e1fs1/e1fs2';
      sys.exit(1);

   if re.match('e1fs1', h10type):
      print '\n<exp>.<dtype> = e1f simulation 1';
      sim=True;
      datadir=os.environ['E1F_SIM2PI_DATADIR1'];
      anadir =os.environ['E1F_SIM2PI_ANADIR1'];
      os.system('cat ${HOME}/ongoing/e1fs1.lst > /tmp/simlst')
   elif re.match('e1fs2', h10type):
      print'\n<exp>.<dtype> = e1f simulation 2';
      sim=True;
      datadir=os.environ['E1F_SIM2PI_DATADIR2'];
      anadir =os.environ['E1F_SIM2PI_ANADIR2'];
      os.system('cat ${HOME}/ongoing/e1fs2.lst > /tmp/simlst')
   elif re.match('e1fd1', h10type):
      print'\n<exp>.<dtype> = e1f experiment 1';
      exp=True;
      datadir=os.environ['E1F_2PI_DATADIR1'];
      anadir =os.environ['E1F_2PI_ANADIR1'];
      os.system('echo "." > /tmp/simlst');
   elif re.match('e1fd2', h10type):
      print '\n<exp>.<dtype> = e1f experiment 2';
      exp=True;
      datadir=os.environ['E1F_2PI_DATADIR2'];
      anadir =os.environ['E1F_2PI_ANADIR2'];
      os.system('echo "." > /tmp/simlst');
   else:
      print 'Could not recognize <exp>.<dtype>=<e1f><exp/sim>';
      sys.exit(1);

   #determine procorder, skim (where applicable)
   if(sim)&(make_yield):
      procorder = 'gProcOrder_simYield';
   elif (sim)&(make_mcyield):
      procorder = 'gProcOrder_mcYield';
   elif (sim)&(make_h10skim):
      procorder = '\"fillskim\"';
      skim='h10skim';
   elif (sim)&(make_h10skim_2piCuts):
      procorder = '\"eid:efid:qskim:pid:fillskim'+monprocs+'\"';
      skim='h10skim_2piCuts';
   elif (exp)&(make_q2wsel):
      procorder = '\"q2wskim:copyh10'+monprocs+'\"';
   elif (exp)&(make_yield):
      procorder = 'gProcOrder_expYield';
         
    

def main():
   # #following determined from argv
   global make_yield,make_mcyield,make_h10skim,make_h10skim_2piCuts,make_q2wsel
   global h10type;
   global skim;
   global sim,exp;
   global bkdirname;
   global nevts;
   global datadir,anadir,procorder;
   
   #get make_<>, h10type, skim(if passed by user), nevt
   parseopts(sys.argv[1:]);
   
   #set up sim,exp,datadir,anadir,procorder,skim
   init();   

   #create logfile
   logfilename='makelog_%s.txt'%datetime.now().strftime('%Y-%m-%d_%H-%M')
   logdir=os.path.join(anadir,'makelog')
   if not os.path.isdir(logdir):
      os.makedirs(logdir)
   logfile=os.path.join(logdir,logfilename);
   sys.stdout = MyOutput(logfile);

   if make_yield: 
      print '### Going to make yield for ###', h10type;
      if bkdirname:
         bkdirname=('bk__yield__'+bkdirname+'__'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))
   elif make_mcyield: 
      print '### Going to make mcyield for ###', h10type;
      if bkdirname:
         bkdirname=('bk__mcyield__'+bkdirname+'__'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))
   elif make_h10skim: print '### Going to make h10skim for ###', h10type;
   elif make_h10skim_2piCuts: print '### Going to make h10skim_2piCuts for ###', h10type;
   elif make_q2wsel: print '### Going to make h10skim_2piCuts for ###', h10type

   fsimlst = open('/tmp/simlst','r')
   for simnum in fsimlst.readlines():
      simnum=simnum.rstrip("\n")#remove blank lines
      print '\n',simnum
      if (make_h10skim|make_h10skim_2piCuts): #make_<skim>
         h10datadir=os.path.join(datadir,simnum,'h10')
         skimdir   =os.path.join(datadir,simnum,skim)
         if not os.path.isdir(skimdir):
            os.makedirs(skimdir)
            print 'going to make', skimdir
            subprocess.call(["sub.runSelector",h10datadir,skimdir,procorder])
         else:
            print skimdir ,'exists! If need to overwrite, then implement in make.py'
      elif (make_q2wsel): #make_q2wsel
         q2wFulldir=os.path.join(datadir,'Q2W__Full','h10')
         q2wseldir =os.path.join(datadir,simnum,'h10')
         if not os.path.isdir(q2wseldir):
            os.makedirs(q2wseldir)
         subprocess.call(["sub.runSelector",q2wFulldir,q2wseldir,procorder])
      else: #make_yield/make_mcyield
         myt=make_tools.YieldTool(datadir,anadir,simnum,skim,procorder,nevts,bkdirname,logfile)
         
         if make_yield:
            myt.setup('cooked')
         elif make_mcyield:
            myt.setup('mc')

         myt.display();
         #sys.exit("Exiting make.py at debug-point")
         myt.run();

   #Add yields
   print '\nNow going to add yields\n'
   if   re.match('e1fs1', h10type):
      if make_yield:   subprocess.call(['addy_e1fs1',bkdirname,logfile])
      if make_mcyield: subprocess.call(['addmcy_e1fs1',bkdirname,logfile])
   elif re.match('e1fs2', h10type):
      if make_yield:   subprocess.call(['addy_e1fs2',bkdirname,logfile])
      if make_mcyield: subprocess.call(['addmcy_e1fs2',bkdirname,logfile])

   print '\n### make.py is done ###\n'
         
if __name__ == "__main__":
   main()
