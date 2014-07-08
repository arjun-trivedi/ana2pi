import ROOT
from rootpy.io import root_open, DoesNotExist

from collections import OrderedDict

import os,sys

class DispYields:
	def __init__(self,q2w):
		self.Q2W=q2w
		self.FEXP=root_open(os.path.join(os.environ['OBS_DIR'],self.Q2W,'yield_exp.root'))
		self.FSIM=root_open(os.path.join(os.environ['OBS_DIR'],self.Q2W,'yield_sim.root'))
		self.ANADIR=os.path.join(os.environ['OBS_DIR'],self.Q2W)
		if not os.path.exists(self.ANADIR):
			sys.exit("%s does not exist!"%self.ANADIR)

	def plot_obs_1D(self,q2wbin,h_dpp,h_rho,h_dzr):
		c=ROOT.TCanvas()
		c.Divide(3,3)
		c.cd(1)
		h_dpp['EXP'][0].DrawNormalized("",1000)
		h_dpp['SIM'][0].DrawNormalized("sames",1000)
		c.cd(2)
		h_rho['EXP'][0].DrawNormalized("",1000)
		h_rho['SIM'][0].DrawNormalized("sames",1000)
		c.cd(3)
		h_dzr['EXP'][0].DrawNormalized("",1000)
		h_dzr['SIM'][0].DrawNormalized("sames",1000)
		c.cd(4)
		h_dpp['EXP'][1].DrawNormalized("",1000)
		h_dpp['SIM'][1].DrawNormalized("sames",1000)
		c.cd(5)
		h_rho['EXP'][1].DrawNormalized("",1000)
		h_rho['SIM'][1].DrawNormalized("sames",1000)
		c.cd(6)
		h_dzr['EXP'][1].DrawNormalized("",1000)
		h_dzr['SIM'][1].DrawNormalized("sames",1000)
		c.cd(7)
		h_dpp['EXP'][2].DrawNormalized("",1000)
		h_dpp['SIM'][2].DrawNormalized("sames",1000)
		c.cd(8)
		h_rho['EXP'][2].DrawNormalized("",1000)
		h_rho['SIM'][2].DrawNormalized("sames",1000)
		c.cd(9)
		h_dzr['EXP'][2].DrawNormalized("",1000)
		h_dzr['SIM'][2].DrawNormalized("sames",1000)
		c.SaveAs("%s/c%s.png"%(self.ANADIR,q2wbin))

	def disp_1D(self):
		"""
		Walk the ROOT file and extract 1D observable histograms. The ROOT file is arranged in a Tree Structure. The
		Observable histograms (obs-hists) are located as files in the following directory-path(dirpath):
		q2wbin/vst/seq/hists
		"""
		#! First get all q2wbin directories from file
		q2ws=[]
		for path,dirs,files in self.FEXP.walk():
			if path=="":continue #! Avoid root path
			path_arr=path.split("/")
			if len(path_arr)==1:
				q2ws.append(path)
		print q2ws

		#! Now get relevant histograms, plot and save them
		for q2w in q2ws:
			h_dpp,h_rho,h_dzr=OrderedDict(),OrderedDict(),OrderedDict()
			for dtyp in ['EXP','SIM']:
				if dtyp=='EXP':
					f=self.FEXP
					col=ROOT.gROOT.ProcessLine("kBlue")
				if dtyp=='SIM':
					f=self.FSIM
					col=ROOT.gROOT.ProcessLine("kRed")
				h_dpp[dtyp]=[f.Get("%s/%s"%(q2w,"VST1/F/h_M1")),f.Get("%s/%s"%(q2w,"VST1/F/h_THETA")),
							 f.Get("%s/%s"%(q2w,"VST1/F/h_ALPHA"))]
				h_rho[dtyp]=[f.Get("%s/%s"%(q2w,"VST2/F/h_M2")),f.Get("%s/%s"%(q2w,"VST2/F/h_THETA")),
							 f.Get("%s/%s"%(q2w,"VST2/F/h_ALPHA"))]
				h_dzr[dtyp]=[f.Get("%s/%s"%(q2w,"VST3/F/h_M2")),f.Get("%s/%s"%(q2w,"VST3/F/h_THETA")),
							 f.Get("%s/%s"%(q2w,"VST3/F/h_ALPHA"))]
				for hists in [h_dpp[dtyp],h_rho[dtyp],h_dzr[dtyp]]:
					#print hists
					for h in hists:
						#print h
						h.SetLineColor(col)
						h.SetMarkerColor(col)
			self.plot_obs_1D(q2w,h_dpp,h_rho,h_dzr)