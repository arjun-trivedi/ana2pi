import ROOT
from rootpy.io import root_open, DoesNotExist

from collections import OrderedDict

import os,sys

import numpy as np
import matplotlib.pyplot as plt

import atlib as atlib

ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool

# Tools 
thntool=THnTool()

class DispYields:
	def __init__(self,simnum='siml'):
		self.SIM_NUM=simnum
		#self.Q2W=q2w
		# self.FEXP=root_open('$HOME/ongoing/mem_test/exp/new-h8-bng/yield_exp.root')
		# self.FSIM=root_open('$HOME/ongoing/mem_test/sim/new-h8-bng/yield_sim.root')
		# self.OUTDIR='/home/trivedia/ongoing/mem_test/obs'
		self.FEXP=root_open(os.path.join(os.environ['OBSDIR'],self.SIM_NUM,'yield_exp.root'))
		self.FSIM=root_open(os.path.join(os.environ['OBSDIR'],self.SIM_NUM,'yield_sim.root'))
		self.OUTDIR=os.path.join(os.environ['OBSDIR'],self.SIM_NUM)
		if not os.path.exists(self.OUTDIR):
			sys.exit("%s does not exist!"%self.OUTDIR)

		self.VSTS=[1,2,3]

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
		
		c.SaveAs("%s/c1D_%s.png"%(self.OUTDIR,q2wbin))

	def disp_1D(self,dtypl=['EXP','SIM'],seql=['T','C','H','F']):
		"""
		Walk the ROOT file and extract 1D observable histograms. 
		"""
		#! 1. First get all q2wbin directories from file
		q2wbinl=self.get_q2wbinl()
		print q2wbinl

		#! 2. Now get relevant histograms
		hVST1,hVST2,hVST3={},{},{}
		q2wbinl_bad={} 
		for q2wbin in q2wbinl:
			print "Processing",q2wbin
			#! First make sure this bin is "good"
			badbin=False
			reason=''
			for vst in self.VSTS:
				ht_R=self.FEXP.Get("%s/VST%d/R/h_M1"%(q2wbin,vst))
				ht_C=self.FEXP.Get("%s/VST%d/C/h_M1"%(q2wbin,vst))
				if ht_R.Integral()==0:
					reason+="ER=0 for VST%d;"%vst
					badbin=True
				elif ht_C.Integral()==0:
					reason+="SA=0 in ER bins for VST%d;"%vst
					badbin=True
			if badbin:
				q2wbinl_bad[q2wbin]=reason
				print "\"Bad\" bin: %s"%q2wbinl_bad[q2wbin]
				continue
			#! Now that bin is good, get hists
			h_dpp,h_rho,h_dzr=OrderedDict(),OrderedDict(),OrderedDict()
			for dtyp in dtypl:
				for seq in seql:
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

					h_dpp[dtyp,seq]=[f.Get("%s/VST1/%s/h_M1"%(q2wbin,seq)),f.Get("%s/VST1/%s/h_THETA"%(q2wbin,seq)),
								 	f.Get("%s/VST1/%s/h_ALPHA"%(q2wbin,seq))]
					h_rho[dtyp,seq]=[f.Get("%s/VST2/%s/h_M2"%(q2wbin,seq)),f.Get("%s/VST2/%s/h_THETA"%(q2wbin,seq)),
								 	f.Get("%s/VST2/%s/h_ALPHA"%(q2wbin,seq))]
					h_dzr[dtyp,seq]=[f.Get("%s/VST3/%s/h_M2"%(q2wbin,seq)),f.Get("%s/VST3/%s/h_THETA"%(q2wbin,seq)),
								 	f.Get("%s/VST3/%s/h_ALPHA"%(q2wbin,seq))]
					for hists in [h_dpp[dtyp,seq],h_rho[dtyp,seq],h_dzr[dtyp,seq]]:
						#print hists
						for h in hists:
							#print h
							h.SetLineColor(col)
							h.SetMarkerColor(col)
							if dtyp=='SIM' and seq=='T':
								 h.SetMarkerStyle(mst)
			self.plot_obs_1D(q2wbin,h_dpp,h_rho,h_dzr)
		print "Following are the \"bad\" q2w bins:"
		for k in q2wbinl_bad:
			print k,q2wbinl_bad[k]

	def disp_integ_yield(self,seq='F'):
		"""
		Walk the ROOT file and plot y(w;q2). 
		"""
		outdir=os.path.join(self.OUTDIR,"integ_yield")
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		#! 1. First of all get q2bng information
		q2bng=self.get_q2bng()
		#print q2bng['BINW']
		print "Yield(w) will be plotted for the following Q2 bins:"
		print q2bng['BINS']

		#! 2. Get all q2wbins
		q2wbinl=self.get_q2wbinl()
		#print q2wbinl

		#! 3. Now put together y[q2bin]=(w,yield)
		y=[{} for i in range(q2bng['NBINS'])]
		for q2wbin in q2wbinl:
			#! Get iq2bin
			q2bin_le=float(q2wbin.split('_')[0].split('-')[0])
			iq2bin=q2bng['BINS_LE'].index(q2bin_le)
			#! Get w, yield
			w=float(q2wbin.split('_')[1].split('-')[0])
			h5=self.FEXP.Get("%s/VST1/%s/h5"%(q2wbin,seq))
			y[iq2bin][w]=thntool.GetIntegral(h5)
		#! Make sure y[q2bin]=(w,yield) are sorted by w
		oy=[{} for i in range(len(y))]
		for iq2bin in range(len(y)):
			oy[iq2bin]=OrderedDict(sorted(y[iq2bin].items()))
			# print iq2bin
			# print y[iq2bin]
			# print oy[iq2bin] 

		#! 4. Now plot
		fig=plt.figure()
		ax=plt.subplot(111)
		clrs=['red','green','cyan','blue','black','yellow','brown','orange']
		for iq2bin in range(len(oy)):
			lbl='[%.2f-%.2f]'%(q2bng['BINS_LE'][iq2bin],q2bng['BINS_UE'][iq2bin])
			clr=clrs[iq2bin]
			ax.scatter(oy[iq2bin].keys(),oy[iq2bin].values(),label=lbl,c=clr)
		ax.set_ylim(0,600000)
		ax.legend()
		fig.savefig('%s/integ_yield.png'%(outdir))	
		#fig.savefig('test.png')		

	# def get_integ_yield(self):
	# 	"""
	# 	Walk the ROOT file and obtain y(seq,w) 
	# 	"""
	# 	#! First get all q2wbin directories from file
	# 	q2ws=self.get_q2wbinl()
	# 	print q2ws

	# 	#! Obtain y(seq,w) for experiment
	# 	#! Example of techincal implementation:
	# 	#! y = { 'C':{w1:y_C,...,wn:y_C},
	# 	#!       'H':{w1:y_H,...,wn:y_H},
	# 	#!       'F':{w1:y_F,...,wn:y_F},
	# 	#!	   }
	# 	y={}
	# 	f=ROOT.TFile(self.FEXP.GetName())
	# 	for seq in ['C','H','F']:
	# 		tmp={}
	# 		for q2w in q2ws:
	# 			w=float(q2w.split('_')[1].split('-')[0])
	# 			h5=f.Get("%s/VST1/%s/h5"%(q2w,seq))
	# 			tmp[w]=thntool.GetIntegral(h5)
	# 		y[seq]=tmp
	# 	return y	

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
		q2ws=self.get_q2wbinl()
		#print "Processing sim_stats for %s:"%self.Q2W
		print q2ws

		ss={'T':[],'R':[],'A':[],'H':[]}
		f=ROOT.TFile(self.FSIM.GetName())
		for seq in ['T','R','A','H']:
			for q2w in q2ws:
				print "Processing sim_stats for %s"%q2w
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
		q2ws=self.get_q2wbinl()
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


	def get_q2wbinl(self):
		"""
		The ROOT file is arranged in a Tree Structure. The
		Observable histograms (obs-hists) are located as files in the following directory-path(dirpath):
		q2wbin/vst/seq/hists
		"""
		q2ws=[]
		i=0 #! for debugging over limited q2ws
		for path,dirs,files in self.FEXP.walk():
			if path=="":continue #! Avoid root path
			path_arr=path.split("/")
			if len(path_arr)==1:
				q2ws.append(path)
				i+=1
			#if i>50: break
		return q2ws

	def get_q2bng(self):
		return self.get_xbng(x='Q2')
	def get_wbng(self):
		return self.get_xbng(x='W')

	def get_xbng(self,x):
		"""
		Gets x=(Q2 or W) binning information from yield_exp.root
		"""

		#! First get all q2w-bins
		q2wbinl=self.get_q2wbinl()
		
		#! 2. From q2wbinl, identify xbng 
		xbins=[]
		for q2wbin in q2wbinl:
			#print q2wbin
			if(x=='Q2'):
				xbins.append(float(q2wbin.split('_')[0].split('-')[0]))
				xbins.append(float(q2wbin.split('_')[0].split('-')[1]))
			elif(x=='W'):
				xbins.append(float(q2wbin.split('_')[1].split('-')[0]))
				xbins.append(float(q2wbin.split('_')[1].split('-')[1]))
			else:
				sys.exit("DispYields::get_xbng() x=%s not recognized"%x)
		xbins=set(xbins) #! Keep only unique entries
		xbins=list(xbins) #! convert 'set' back to 'list'
		xbins=sorted(xbins) #! Order entries
		#print xbins
		xbins_le=xbins[:-1]
		xbins_ue=xbins[1:]
		#print xbins_le
		#print xbins_ue
		xmin=min(xbins)
		xmax=max(xbins)
		xbinw=xbins[1]-xbins[0]
		nxbins=len(xbins_le)
		xbng={'MIN':xmin,'MAX':xmax,'BINW':xbinw,'NBINS':nxbins,
		   'BINS_LE':xbins_le,'BINS_UE':xbins_ue,'BINS':xbins}
		return xbng