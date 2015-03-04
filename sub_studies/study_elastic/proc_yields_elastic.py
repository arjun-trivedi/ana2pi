import os
from collections import OrderedDict
import ROOT

H2_DIM=OrderedDict([('THETA',0),('PHI',1)])
def run():
        #! Get all input delast files
	FIN={}
	FIN['ER']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_EXP'],'delastR.root'))
	FIN['SR']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml','delastR.root'))
	FIN['ST']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml','delastT.root'))
	print FIN['ER'].GetName()
	
	FOUT=ROOT.TFile(os.path.join(os.environ['OBSDIR_ELASTIC'],'yield.root'),"RECREATE")

	#! Begin calculation
	hN2={}
	for seq in ['ST','SR','SA','ER','EC']:
		print seq
		if seq=='ST' or seq=='SR' or seq=='ER':
			print "passed"
			fin=FIN[seq]
			hN2[seq]=fin.Get("delast/yield")
			FOUT.mkdir(seq).cd()
			hN2[seq].Write()
			h2=hN2[seq].Projection(H2_DIM['PHI'],H2_DIM['THETA'])
			h2.SetName("h2")
			h2.Write()
			
