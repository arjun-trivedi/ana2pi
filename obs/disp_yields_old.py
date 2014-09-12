import ROOT
from rootpy.io import root_open, DoesNotExist

from collections import OrderedDict

import os,sys

import numpy as np

ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool

# Tools 
thntool=THnTool()

class DispYields:
	def __init__(self,q2w,sim_num='siml'):
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
		# h_dpp['EXP','F'][0].Draw()
		# h_dpp['EXP','H'][0].Draw("sames")
		hexp=h_dpp['EXP','F'][0].DrawNormalized("",1000)
		hsim=h_dpp['SIM','F'][0].DrawNormalized("sames",1000)
		hexp.SetMinimum(0.)
		hsim.SetMinimum(0.)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		hF=hexp.Clone()
		hF.Divide(h_dpp['EXP','F'][0])
		h_dpp['EXP','C'][0].Multiply(hF)
		h_dpp['EXP','C'][0].Draw("sames")
		h_dpp['EXP','H'][0].Multiply(hF)
		h_dpp['EXP','H'][0].Draw("sames")
		#h_dpp['EXP','H'][0].DrawNormalized("sames",1000)
		h_dpp['SIM','T'][0].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(2)
		# h_rho['EXP','F'][0].Draw()
		# h_rho['EXP','H'][0].Draw("sames")
		hexp=h_rho['EXP','F'][0].DrawNormalized("",1000)
		hsim=h_rho['SIM','F'][0].DrawNormalized("sames",1000)
		hexp.SetMinimum(0.)
		hsim.SetMinimum(0.)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		hF=hexp.Clone()
		hF.Divide(h_rho['EXP','F'][0])
		h_rho['EXP','C'][0].Multiply(hF)
		h_rho['EXP','C'][0].Draw("sames")
		h_rho['EXP','H'][0].Multiply(hF)
		h_rho['EXP','H'][0].Draw("sames")
		#h_rho['EXP','H'][0].DrawNormalized("sames",1000)
		h_rho['SIM','T'][0].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(3)
		# h_dzr['EXP','F'][0].Draw()
		# h_dzr['EXP','H'][0].Draw("sames")
		hexp=h_dzr['EXP','F'][0].DrawNormalized("",1000)
		hsim=h_dzr['SIM','F'][0].DrawNormalized("sames",1000)
		hexp.SetMinimum(0.)
		hsim.SetMinimum(0.)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		hF=hexp.Clone()
		hF.Divide(h_dzr['EXP','F'][0])
		h_dzr['EXP','C'][0].Multiply(hF)
		h_dzr['EXP','C'][0].Draw("sames")
		h_dzr['EXP','H'][0].Multiply(hF)
		h_dzr['EXP','H'][0].Draw("sames")
		#h_dzr['EXP','H'][0].DrawNormalized("sames",1000)
		h_dzr['SIM','T'][0].DrawNormalized("sames",1000)
		# pad.Update()

		pad=c.cd(4)
		# h_dpp['EXP','F'][1].Draw()
		# h_dpp['EXP','H'][1].Draw("sames")
		hexp=h_dpp['EXP','F'][1].DrawNormalized("",1000)
		hsim=h_dpp['SIM','F'][1].DrawNormalized("sames",1000)
		hexp.SetMinimum(0.)
		hsim.SetMinimum(0.)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		hF=hexp.Clone()
		hF.Divide(h_dpp['EXP','F'][1])
		h_dpp['EXP','C'][1].Multiply(hF)
		h_dpp['EXP','C'][1].Draw("sames")
		h_dpp['EXP','H'][1].Multiply(hF)
		h_dpp['EXP','H'][1].Draw("sames")
		#h_dpp['EXP','H'][1].DrawNormalized("sames",1000)
		h_dpp['SIM','T'][1].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(5)
		# h_rho['EXP','F'][1].Draw()
		# h_rho['EXP','H'][1].Draw("sames")
		hexp=h_rho['EXP','F'][1].DrawNormalized("",1000)
		hsim=h_rho['SIM','F'][1].DrawNormalized("sames",1000)
		hexp.SetMinimum(0.)
		hsim.SetMinimum(0.)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		hF=hexp.Clone()
		hF.Divide(h_rho['EXP','F'][1])
		h_rho['EXP','C'][1].Multiply(hF)
		h_rho['EXP','C'][1].Draw("sames")
		h_rho['EXP','H'][1].Multiply(hF)
		h_rho['EXP','H'][1].Draw("sames")
		#h_rho['EXP','H'][1].DrawNormalized("sames",1000)
		h_rho['SIM','T'][1].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(6)
		# h_dzr['EXP','F'][1].Draw()
		# h_dzr['EXP','H'][1].Draw("sames")
		hexp=h_dzr['EXP','F'][1].DrawNormalized("",1000)
		hsim=h_dzr['SIM','F'][1].DrawNormalized("sames",1000)
		hexp.SetMinimum(0.)
		hsim.SetMinimum(0.)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		hF=hexp.Clone()
		hF.Divide(h_dzr['EXP','F'][1])
		h_dzr['EXP','C'][1].Multiply(hF)
		h_dzr['EXP','C'][1].Draw("sames")
		h_dzr['EXP','H'][1].Multiply(hF)
		h_dzr['EXP','H'][1].Draw("sames")
		#h_dzr['EXP','H'][1].DrawNormalized("sames",1000)
		h_dzr['SIM','T'][1].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(7)
		# h_dpp['EXP','F'][2].Draw()
		# h_dpp['EXP','H'][2].Draw("sames")
		hexp=h_dpp['EXP','F'][2].DrawNormalized("",1000)
		hsim=h_dpp['SIM','F'][2].DrawNormalized("sames",1000)
		hexp.SetMinimum(0.)
		hsim.SetMinimum(0.)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		hF=hexp.Clone()
		hF.Divide(h_dpp['EXP','F'][2])
		h_dpp['EXP','C'][2].Multiply(hF)
		h_dpp['EXP','C'][2].Draw("sames")
		h_dpp['EXP','H'][2].Multiply(hF)
		h_dpp['EXP','H'][2].Draw("sames")
		#h_dpp['EXP','H'][2].DrawNormalized("sames",1000)
		h_dpp['SIM','T'][2].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(8)
		# h_rho['EXP','F'][2].Draw()
		# h_rho['EXP','H'][2].Draw("sames")
		hexp=h_rho['EXP','F'][2].DrawNormalized("",1000)
		hsim=h_rho['SIM','F'][2].DrawNormalized("sames",1000)
		hexp.SetMinimum(0.)
		hsim.SetMinimum(0.)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		hF=hexp.Clone()
		hF.Divide(h_rho['EXP','F'][2])
		h_rho['EXP','C'][2].Multiply(hF)
		h_rho['EXP','C'][2].Draw("sames")
		h_rho['EXP','H'][2].Multiply(hF)
		h_rho['EXP','H'][2].Draw("sames")
		#h_rho['EXP','H'][2].DrawNormalized("sames",1000)
		h_rho['SIM','T'][2].DrawNormalized("sames",1000)
		pad.Update()

		pad=c.cd(9)
		# h_dzr['EXP','F'][2].Draw()
		# h_dzr['EXP','H'][2].Draw("sames")
		hexp=h_dzr['EXP','F'][2].DrawNormalized("",1000)
		hsim=h_dzr['SIM','F'][2].DrawNormalized("sames",1000)
		hexp.SetMinimum(0.)
		hsim.SetMinimum(0.)
		maximum=hexp.GetMaximum()
		if hsim.GetMaximum()>hexp.GetMaximum():
			maximum=hsim.GetMaximum()
		hexp.SetMaximum(maximum+10)
		hsim.SetMaximum(maximum+10)
		hF=hexp.Clone()
		hF.Divide(h_dzr['EXP','F'][2])
		h_dzr['EXP','C'][2].Multiply(hF)
		h_dzr['EXP','C'][2].Draw("sames")
		h_dzr['EXP','H'][2].Multiply(hF)
		h_dzr['EXP','H'][2].Draw("sames")
		#h_dzr['EXP','H'][2].DrawNormalized("sames",1000)
		h_dzr['SIM','T'][2].DrawNormalized("sames",1000)
		pad.Update()
		
		c.SaveAs("%s/c1D_%s.png"%(self.ANADIR,q2wbin))

	def disp_1D(self):
		"""
		Walk the ROOT file and extract 1D observable histograms. 
		"""
		#! First get all q2wbin directories from file
		q2ws=self.get_q2ws()
		print q2ws

		#! Now get relevant histograms, plot and save them
		for q2w in q2ws:
			h_dpp,h_rho,h_dzr=OrderedDict(),OrderedDict(),OrderedDict()
			for dtyp in ['EXP','SIM']:
				for seq in ['T','C','H','F']:
					if dtyp=='EXP' and seq=='T': continue
					if dtyp=='SIM' and (seq=='H' or seq=='C'): continue

					if dtyp=='EXP':
						f=self.FEXP
						if seq=='C':
							col=ROOT.gROOT.ProcessLine("kCyan")
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

	def get_integ_yield(self):
		"""
		Walk the ROOT file and obtain y(seq,w) 
		"""
		#! First get all q2wbin directories from file
		q2ws=self.get_q2ws()
		print q2ws

		#! Obtain y(seq,w) for experiment
		#! Example of techincal implementation:
		#! y = { 'C':{w1:y_C,...,wn:y_C},
		#!       'H':{w1:y_H,...,wn:y_H},
		#!       'F':{w1:y_F,...,wn:y_F},
		#!	   }
		y={}
		f=ROOT.TFile(self.FEXP.GetName())
		for seq in ['C','H','F']:
			tmp={}
			for q2w in q2ws:
				w=float(q2w.split('_')[1].split('-')[0])
				h5=f.Get("%s/VST1/%s/h5"%(q2w,seq))
				tmp[w]=thntool.GetIntegral(h5)
			y[seq]=tmp
		return y	

	def get_sim_stats(self):
		"""
		Walk the ROOT file and obtain simstats(ss) for a h5 in a Q2-W bin:
		ss={'T':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'R':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'A':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'H':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]]}

		where:
			+ nbins=number of filled bins in a h5
			+ N=sum({n_i}) were n_i=events per bin (number of events in a h5)
			+ mu=average({n_i} (average of number of events per bin)
			+ sg=(RMS({n_i}) (RMS of number of events per bin)
		"""
		#! First get all q2wbin directories from file
		q2ws=self.get_q2ws()
		print "Processing sim_stats for %s:"%self.Q2W
		print q2ws

		ss={'T':[],'R':[],'A':[],'H':[]}
		f=ROOT.TFile(self.FSIM.GetName())
		for seq in ['T','R','A','H']:
			for q2w in q2ws:
				#! Determine q2,w
				q2bin=q2w.split('_')[0]
				wbin=q2w.split('_')[1]
				#print q2bin,wbin
				q2=float(q2bin.split('-')[0])
				w=float(wbin.split('-')[0])
				#print q2,w
				#! Determine nbins,N,mu,sg for this q2,w
				h5=f.Get("%s/VST1/%s/h5"%(q2w,seq))
				nbins=thntool.GetNbinsNotEq0(h5)
				N=thntool.GetIntegral(h5)
				binc_stats=np.zeros(2,'f')
				thntool.GetBinContentDistStats(h5,binc_stats)
				mu=binc_stats[0]
				sg=binc_stats[1]
				ss[seq].append([q2,w,nbins,N,mu,sg])
			# #! Compute average
			# ss[seq].append(nevts/len(q2ws))
			# ss[seq].append(nbins/len(q2ws))
		return ss

	def get_sim_stats_commonbins(self):
		"""
		The same as get_sim_stats(), but with the difference that mu,sg
		are obtained only for R-PS bins. I am currently not using this since
		this process takes impractically long.

		Walk the ROOT file and obtain simstats(ss) for a h5 in a Q2-W bin:
		ss={'T':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'R':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'A':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'H':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]]}

		where:
			+ nbins=number of filled bins in a h5
			+ N=sum({n_i}) were n_i=events per bin (number of events in a h5)
			+ mu=average({n_i}) (average of number of events per bin)
				+ Note, average computed over only R-PS bins
			+ sg=(RMS({n_i}) (RMS of number of events per bin)
				+ Note, average computed over only R-PS bins
		"""
		#! First get all q2wbin directories from file
		q2ws=self.get_q2ws()
		print "Processing sim_stats for %s:"%self.Q2W
		print q2ws

		ss={'T':[],'R':[],'A':[],'H':[]}
		f=ROOT.TFile(self.FSIM.GetName())
		for q2w in q2ws:
			print "Processing %s..."%q2w
			#! Determine q2,w
			q2bin=q2w.split('_')[0]
			wbin=q2w.split('_')[1]
			#print q2bin,wbin
			q2=float(q2bin.split('-')[0])
			w=float(wbin.split('-')[0])
			#print q2,w
			#! First get all h5s
			h5={}
			for seq in ['T','R','A','H']:
				h5[seq]=f.Get("%s/VST1/%s/h5"%(q2w,seq))
			#! Now get simstats	
			for seq in ['T','R','A','H']:
				#! Determine nbins,N,mu,sg for this q2,w
				nbins=thntool.GetNbinsNotEq0(h5[seq])
				N=thntool.GetIntegral(h5[seq])
				binc_stats=np.zeros(2,'f')
				print "here"
				thntool.GetBinContentDistStatsCommonBins(h5[seq],h5['R'],binc_stats)
				print "here1"
				mu=binc_stats[0]
				sg=binc_stats[1]
				ss[seq].append([q2,w,nbins,N,mu,sg])
		return ss


	def get_q2ws(self):
		"""
		The ROOT file is arranged in a Tree Structure. The
		Observable histograms (obs-hists) are located as files in the following directory-path(dirpath):
		q2wbin/vst/seq/hists
		"""
		q2ws=[]
		for path,dirs,files in self.FEXP.walk():
			if path=="":continue #! Avoid root path
			path_arr=path.split("/")
			if len(path_arr)==1:
				q2ws.append(path)
		return q2ws