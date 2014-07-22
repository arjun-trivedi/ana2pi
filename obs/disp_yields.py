import ROOT
from rootpy.io import root_open, DoesNotExist

from collections import OrderedDict

import os,sys

class DispYields:
	def __init__(self,q2w,sim_num='sim1'):
		self.SIM_NUM=sim_num
		self.Q2W=q2w
		self.FEXP=root_open(os.path.join(os.environ['OBS_DIR'],self.SIM_NUM,self.Q2W,'yield_exp.root'))
		self.FSIM=root_open(os.path.join(os.environ['OBS_DIR'],self.SIM_NUM,self.Q2W,'yield_sim.root'))
		self.ANADIR=os.path.join(os.environ['OBS_DIR'],self.SIM_NUM,self.Q2W)
		if not os.path.exists(self.ANADIR):
			sys.exit("%s does not exist!"%self.ANADIR)

	def plot_obs_1D(self,q2wbin,h_dpp,h_rho,h_dzr):
		c=ROOT.TCanvas()
		c.Divide(3,3)
		pad=c.cd(1)
		hexp=h_dpp['EXP','F'][0].DrawNormalized("",1000)
		hsim=h_dpp['SIM','F'][0].DrawNormalized("sames",1000)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		h_dpp['EXP','H'][0].DrawNormalized("sames",1000)
		h_dpp['SIM','T'][0].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(2)
		hexp=h_rho['EXP','F'][0].DrawNormalized("",1000)
		hsim=h_rho['SIM','F'][0].DrawNormalized("sames",1000)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		h_rho['EXP','H'][0].DrawNormalized("sames",1000)
		h_rho['SIM','T'][0].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(3)
		hexp=h_dzr['EXP','F'][0].DrawNormalized("",1000)
		hsim=h_dzr['SIM','F'][0].DrawNormalized("sames",1000)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		h_dzr['EXP','H'][0].DrawNormalized("sames",1000)
		h_dzr['SIM','T'][0].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(4)
		hexp=h_dpp['EXP','F'][1].DrawNormalized("",1000)
		hsim=h_dpp['SIM','F'][1].DrawNormalized("sames",1000)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		h_dpp['EXP','H'][1].DrawNormalized("sames",1000)
		h_dpp['SIM','T'][1].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(5)
		hexp=h_rho['EXP','F'][1].DrawNormalized("",1000)
		hsim=h_rho['SIM','F'][1].DrawNormalized("sames",1000)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		h_rho['EXP','H'][1].DrawNormalized("sames",1000)
		h_rho['SIM','T'][1].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(6)
		hexp=h_dzr['EXP','F'][1].DrawNormalized("",1000)
		hsim=h_dzr['SIM','F'][1].DrawNormalized("sames",1000)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		h_dzr['EXP','H'][1].DrawNormalized("sames",1000)
		h_dzr['SIM','T'][1].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(7)
		hexp=h_dpp['EXP','F'][2].DrawNormalized("",1000)
		hsim=h_dpp['SIM','F'][2].DrawNormalized("sames",1000)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		h_dpp['EXP','H'][2].DrawNormalized("sames",1000)
		h_dpp['SIM','T'][2].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(8)
		hexp=h_rho['EXP','F'][2].DrawNormalized("",1000)
		hsim=h_rho['SIM','F'][2].DrawNormalized("sames",1000)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		h_rho['EXP','H'][2].DrawNormalized("sames",1000)
		h_rho['SIM','T'][2].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(9)
		hexp=h_dzr['EXP','F'][2].DrawNormalized("",1000)
		hsim=h_dzr['SIM','F'][2].DrawNormalized("sames",1000)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		h_dzr['EXP','H'][2].DrawNormalized("sames",1000)
		h_dzr['SIM','T'][2].DrawNormalized("sames",1000)
		pad.Update()
		
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
				for seq in ['T','H','F']:
					if dtyp=='EXP' and seq=='T': continue
					if dtyp=='SIM' and seq=='H': continue

					if dtyp=='EXP':
						f=self.FEXP
						if seq=='F':
							col=ROOT.gROOT.ProcessLine("kBlue")
						if seq=='H':
							col=ROOT.gROOT.ProcessLine("kBlack")
					if dtyp=='SIM':
						f=self.FSIM
						if seq=='F':
							col=ROOT.gROOT.ProcessLine("kRed")
						if seq=='T':
							col=ROOT.gROOT.ProcessLine("kGreen")
							mst=ROOT.gROOT.ProcessLine("kPlus")

					h_dpp[dtyp,seq]=[f.Get("%s/VST1/%s/h_M1"%(q2w,seq)),f.Get("%s/VST1/%s/h_THETA"%(q2w,seq)),
								 	f.Get("%s/VST1/%s/h_ALPHA"%(q2w,seq))]
					h_rho[dtyp,seq]=[f.Get("%s/VST2/%s/h_M2"%(q2w,seq)),f.Get("%s/VST2/%s/h_THETA"%(q2w,seq)),
								 	f.Get("%s/VST2/%s/h_ALPHA"%(q2w,seq))]
					h_dzr[dtyp,seq]=[f.Get("%s/VST3/%s/h_M2"%(q2w,seq)),f.Get("%s/VST3/%s/h_THETA"%(q2w,seq)),
								 	f.Get("%s/VST3/%s/h_ALPHA"%(q2w,seq))]
					for hists in [h_dpp[dtyp,seq],h_rho[dtyp,seq],h_dzr[dtyp,seq]]:
						#print hists
						for h in hists:
							#print h
							h.SetLineColor(col)
							h.SetMarkerColor(col)
							if dtyp=='SIM' and seq=='T':
								 h.SetMarkerStyle(mst)
			self.plot_obs_1D(q2w,h_dpp,h_rho,h_dzr)