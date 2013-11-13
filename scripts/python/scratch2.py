def plotR2_D(h5):
	norm = 50000
	hR2 = {}
	hR2[(THETA,POS,D)] = h5[(EC,POS,D)].Projection(THETA)
	hR2[(THETA,POS,D)].Scale(1/math.pi)
	hR2[(THETA,POS,D)].Scale(1/norm)
	hR2[(THETA,NEG,D)] = h5[(EC,NEG,D)].Projection(THETA)
	hR2[(THETA,NEG,D)].Scale(1/math.pi)
	hR2[(THETA,NEG,D)].Scale(1/norm)
	hR2[(THETA,AVG,D)] = hR2[(THETA,POS,D)].Clone("avg")
	hR2[(THETA,AVG,D)].Add(hR2[(THETA,NEG,D)])
	hR2[(THETA,AVG,D)].Scale(0.5)
	hR2[(THETA,AVG,D)].SetMinimum(-0.003)
	hR2[(THETA,AVG,D)].SetMaximum(0.003)
	hR2[(THETA,AVG,D)].SetLineColor(gROOT.ProcessLine("kMagenta"));
	hR2[(THETA,AVG,D)].SetMarkerStyle(gROOT.ProcessLine("kFullCircle"));
	#Make Titles nice
	hR2[(THETA,AVG,D)].SetTitle("")
	pt = TPaveText(0.3, 0.85, 0.7, 1.0, "NDC")
	q2wt = pt.AddText('[Q^{2}][W] = %s'%q2wdir.GetName())
	q2wt.SetTextColor(gROOT.ProcessLine("kBlue"))
	vart = pt.AddText(("D^{%s} vs. %s")%(vartitle[0][THETA],vartitle[0][THETA]));
	vart.SetTextSize(0.05);

	
	cR2 = {}
	l = TLine(0,0,180,0)
	cR2[(THETA,AVG, D)] = TCanvas("RvVar", "RvVar")
	hR2[(THETA,AVG, D)].Draw("ep")
	l.Draw("same")
	pt.Draw()
	cR2[(THETA,AVG, D)].SaveAs( ('%s/%s.png')%(outdir,cR2[(THETA,AVG, D)].GetName()))