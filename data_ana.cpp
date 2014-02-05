#include "data_ana.h"
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TMath.h>
#include "particle_constants.h"

#include <iostream>
using namespace std;

using namespace TMath;
using namespace AnalysisConstants;

DataAna::DataAna()
{
	
}

DataAna::~DataAna()
{
}

void DataAna::Clear()
{
	opart = h10idxE = h10idxP = h10idxPip = h10idxPim = -1;
	top = 0;
	eid.Clear();
	efid.Clear();
	skimq.Clear();
	mom.Clear();
	pid.Clear();
	eKin.Clear();
	eKin_mc.Clear();
	d2pi.Clear();
	d2pi_mc.Clear();
	dElast.Clear();
	dElast_ST.Clear();
}

void DataAna::makeHistsEid(TObjArray** hists, TDirectory* dirout)
{
	//TObjArray *ret = new TObjArray(5);
	//TObjArray* ret[NSECTORS];
	//TDirectory* direid = dirout->mkdir("eid");
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
		
		hists[iSector] = new TObjArray(7);
		
		hists[iSector]->Add(new TH2F("heoutVein",TString::Format("#DeltaE_{out} vs. #DeltaE_{in}(sector%d)", iSector),150, 0, 1, 150, 0, 1));
		hists[iSector]->Add(new TH2F("hsfoutVsfin", TString::Format("SF_{out} vs. SF_{in}(sector%d)", iSector),100, 0, 0.5, 100, 0, 0.5));
		hists[iSector]->Add(new TH2F("hsftotVp", TString::Format("SF_{tot} vs. p(sector%d)", iSector), 160, 0, 5, 100, 0, 0.5));
		hists[iSector]->Add(new TH1F("hnphe", TString::Format("nphe(sector%d)", iSector), 100, 0, 300));
		hists[iSector]->Add(new TH2F("hCCthetaVCCseg" ,TString::Format("CC_{#theta} vs. CC_{seg}(sector%d)", iSector), 20, 0, 20, 100, 0, 50));
		hists[iSector]->Add(new TH2F("hbetaVp",TString::Format("#beta vs p(sector%d)", iSector),100, 0, 4, 100, 0, 1.5));
		hists[iSector]->Add(new TH2F("hdtIfVp",TString::Format("#Deltat vs p(sector%d)", iSector),100, 0, 4, 100, -80, 120));
	}
	/*for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		Info("DataAna::makeHistsEid()", "sector%d = %d", iSector, hists[iSector]->GetEntriesFast());
		//Info("DataAna::makeHistsEid()", "sector%d = %s", iSector, hists[iSector]->At(0)->GetTitle());
		Info("DataAna::makeHistsEid()", "sector%d = %s", iSector, hists[iSector]->At(0)->GetTitle());
		Info("DataAna::makeHistsEid()", "sector%d = %s", iSector, hists[iSector]->At(1)->GetTitle());
		Info("DataAna::makeHistsEid()", "sector%d = %s", iSector, hists[iSector]->At(2)->GetTitle());
		Info("DataAna::makeHistsEid()", "sector%d = %s", iSector, hists[iSector]->At(3)->GetTitle());
		Info("DataAna::makeHistsEid()", "sector%d = %s", iSector, hists[iSector]->At(4)->GetTitle());
		Info("DataAna::makeHistsEid()", "sector%d = %s", iSector, hists[iSector]->At(5)->GetTitle());
		Info("DataAna::makeHistsEid()", "sector%d = %s", iSector, hists[iSector]->At(6)->GetTitle());
	}*/
	//return ret;
	return;
}

void DataAna::makeHistsEFid(TObjArray** hists, TDirectory* dirout)
{
	//TObjArray *ret = new TObjArray(6);
	//TObjArray* ret[7];
	//TDirectory* direfid = dirout->mkdir("efid");
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(1);
		
		if (iSector==0)	hists[iSector]->Add(new TH2F("hephiVetheta","#phi(sector 0) vs. #theta", 100, 0, 60, 100, -30, 360));
		if (iSector==1) hists[iSector]->Add(new TH2F("hephiVetheta","#phi(sector 1) vs. #theta", 100, 0, 60, 100, -30, 30));
		if (iSector==2) hists[iSector]->Add(new TH2F("hephiVetheta","#phi(sector 2) vs. #theta", 100, 0, 60, 100, 30, 90));
		if (iSector==3)	hists[iSector]->Add(new TH2F("hephiVetheta","#phi(sector 3) vs. #theta", 100, 0, 60, 100, 90, 150));
		if (iSector==4) hists[iSector]->Add(new TH2F("hephiVetheta","#phi(sector 4) vs. #theta", 100, 0, 60, 100, 150, 210));
		if (iSector==5) hists[iSector]->Add(new TH2F("hephiVetheta","#phi(sector 5) vs. #theta", 100, 0, 60, 100, 210, 270));
		if (iSector==6) hists[iSector]->Add(new TH2F("hephiVetheta","#phi(sector 6) vs. #theta", 100, 0, 60, 100, 270, 330));
		hists[iSector]->Add(new TH2F("hescyVescx", "SCx vs. SCy", 200,0,500,200,-250,250));
		hists[iSector]->Add(new TH2F("heecyVeecx", "ECx vs. ECy", 200,-500,500,200,-500,500));
	}
	//return ret;
	return;
}

void DataAna::makeHistsMomCor(TObjArray** hists, TDirectory* dirout)
{
	//TObjArray *ret = new TObjArray(5);
	//TObjArray* ret[7];
	//TDirectory* dirmomcor = dirout->mkdir("momcor");
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(5);
		
		hists[iSector]->Add(new TH2F("hdpVp",TString::Format("#Deltap vs. p(sector%d)", iSector),550,0,5.5,160,-0.08,0.08));
		hists[iSector]->Add(new TH1F("hdcx", TString::Format("#Deltacx(sector%d)", iSector),60,-0.01,0.01));
		hists[iSector]->Add(new TH1F("hdcy", TString::Format("#Deltacy(sector%d)", iSector), 60,-0.01,0.01));
		hists[iSector]->Add(new TH1F("hdcz", TString::Format("#Deltacz(sector%d)(sector%d)", iSector), 60,-0.01,0.01));
		hists[iSector]->Add(new TH1F("hdp", TString::Format("#Deltap(sector%d)(sector%d)", iSector), 160,-0.08,0.08));
	}
	//return ret;
	return;
}

void DataAna::makeHistsPid(TObjArray** hists, TDirectory* dirout)
{
	//TObjArray *ret = new TObjArray(21);
	TObjArray* ret[7];
	//TDirectory* dirpid = dirout->mkdir("pid");
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(21);
	
		hists[iSector]->Add(new TH2F("hP_betaVp",TString::Format("#beta vs p for idtfd proton(sector%d)", iSector),100, 0, 4, 100, 0, 1.5));
		hists[iSector]->Add(new TH2F("hP_dtVp",TString::Format("#Deltat vs p for idtfd proton(sector%d)", iSector),100, 0, 4, 200, -40, 120));
		hists[iSector]->Add(new TH2F("hP_eoutVein",TString::Format("#DeltaE_{out} vs. #DeltaE_{in} for idtfd p(sector%d)", iSector),150, 0, 1, 150, 0, 1));
		hists[iSector]->Add(new TH2F("hP_sfoutVsfin",TString::Format("SF_{out} vs. SF_{in} for idtfd p(sector%d)", iSector),100, 0, 0.5, 100, 0, 0.5));
		hists[iSector]->Add(new TH2F("hP_sftotVp",TString::Format("SF_{tot} vs. p for idtfd p(sector%d)", iSector), 160, 0, 5, 100, 0, 0.5));
		hists[iSector]->Add(new TH1F("hP_nphe", TString::Format("nphe for idtfd p(sector%d)", iSector), 100, 0, 300));
		hists[iSector]->Add(new TH2F("hP_CCthetaVCCseg",TString::Format("CC_{#theta} vs. CC_{seg} for idtf p(sector%d)", iSector), 20, 0, 20, 100, 0, 50));
		
		hists[iSector]->Add(new TH2F("hPip_betaVp", TString::Format("#beta vs p for idtfd #pi^{+}(sector%d)", iSector),100, 0, 4, 100, 0, 1.5));
		hists[iSector]->Add(new TH2F("hPip_dtVp", TString::Format("#Deltat vs p for idtfd #pi^{+}(sector%d)", iSector),100, 0, 4, 200, -40, 120));
		hists[iSector]->Add(new TH2F("hPip_eoutVein", TString::Format("#DeltaE_{out} vs. #DeltaE_{in} for idtfd #pi^{+}(sector%d)", iSector),150, 0, 1, 150, 0, 1));
		hists[iSector]->Add(new TH2F("hPip_sfoutVsfin", TString::Format("SF_{out} vs. SF_{in} for idtfd #pi^{+}(sector%d)", iSector),100, 0, 0.5, 100, 0, 0.5));
		hists[iSector]->Add(new TH2F("hPip_sftotVp", TString::Format("SF_{tot} vs. p for idtfd #pi^{+}(sector%d)", iSector), 160, 0, 5, 100, 0, 0.5));
		hists[iSector]->Add(new TH1F("hPip_nphe", TString::Format("nphe for idtfd #pi^{+}(sector%d)", iSector), 100, 0, 300));
		hists[iSector]->Add(new TH2F("hPip_CCthetaVCCseg", TString::Format("CC_{#theta} vs. CC_{seg} for idtf #pi^{+}(sector%d)", iSector), 20, 0, 20, 100, 0, 50));
		
		hists[iSector]->Add(new TH2F("hPim_betaVp", TString::Format("#beta vs p for idtfd #pi^{-}(sector%d)", iSector),100, 0, 4, 100, 0, 1.5));
		hists[iSector]->Add(new TH2F("hPim_dtVp", TString::Format("#Deltat vs p for idtfd #pi^{-}(sector%d)", iSector),100, 0, 4, 200, -40, 120));
		hists[iSector]->Add(new TH2F("hPim_eoutVein", TString::Format("#DeltaE_{out} vs. #DeltaE_{in} for idtfd #pi^{-}(sector%d)", iSector),150, 0, 1, 150, 0, 1));
		hists[iSector]->Add(new TH2F("hPim_sfoutVsfin", TString::Format("SF_{out} vs. SF_{in} for idtfd #pi^{-}(sector%d)", iSector),100, 0, 0.5, 100, 0, 0.5));
		hists[iSector]->Add(new TH2F("hPim_sftotVp", TString::Format("SF_{tot} vs. p for idtfd #pi^{-}(sector%d)", iSector), 160, 0, 5, 100, 0, 0.5));
		hists[iSector]->Add(new TH1F("hPim_nphe", TString::Format("nphe for idtfd #pi^{-}(sector%d)", iSector), 100, 0, 300));
		hists[iSector]->Add(new TH2F("hPim_CCthetaVCCseg", TString::Format("CC_{#theta} vs. CC_{seg} for idtf #pi^{-}(sector%d)", iSector), 20, 0, 20, 100, 0, 50));
	}
	//return ret;
	return;
}

void DataAna::makeHistsEkin(TObjArray** hists, TDirectory* dirout)
{
	//TObjArray *ret = new TObjArray(4);
	//TObjArray* ret[NSECTORS];
	//TDirectory* direkin = dirout->mkdir("ekin");
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(6);
		
		hists[iSector]->Add(new TH1F("helastic",TString::Format("W_{elastic}(sector%d)", iSector),100,0.7,1.2));
		hists[iSector]->Add(new TH2F("hq2Vw",TString::Format("Q^{2} vs. W(sector%d)", iSector),600,0,3.6,100,0,6)); //to match Evgeny;s q2VW bng (Q2/bin, W/bin) = (0.06[GeV^2/bin], 0.006[GeV/bin])
		hists[iSector]->Add(new TH1F("hphiq",TString::Format("#phi^{*}(sector%d)", iSector),360,-30,330));
		hists[iSector]->Add(new TH2F("he1Vtheta1",TString::Format("E\' vs. #theta\'(sector%d)", iSector),120,0,60,225,0,5.5));
		hists[iSector]->Add(new TH2F("hq2Vtheta1",TString::Format("Q^{2} vs. #theta\'(sector%d)", iSector),120,0,60,100,0,6));
		hists[iSector]->Add(new TH2F("hq2Ve1",TString::Format("Q^{2} vs. E\'(sector%d)", iSector),225,0,5.5,100,0,6));
	}
	//return ret;
	return;
}

TObjArray* DataAna::makeHistsEkin()
{
	TObjArray *ret = new TObjArray(6);
	ret->Add(new TH1F("helastic","W_{elastic}",100,0.7,1.2));
	ret->Add(new TH2F("hq2Vw","Q^{2} vs. W",600,0,3.6,100,0,6)); //to match Evgeny;s q2VW bng (Q2/bin, W/bin) = (0.06[GeV^2/bin], 0.006[GeV/bin])
	ret->Add(new TH1F("hphiq","#phi^{*}",360,-30,330));
	ret->Add(new TH2F("he1Vtheta1","E\' vs. #theta\'",120,0,60,225,0,5.5));
	ret->Add(new TH2F("hq2Vtheta1","Q^{2} vs. #theta\'",120,0,60,100,0,6));
	ret->Add(new TH2F("hq2Ve1","Q^{2} vs. E\'",225,0,5.5,100,0,6));
	return ret;
}

TObjArray* DataAna::makeHistsMM()
{
	TObjArray *ret = new TObjArray(4);
	ret->Add(new TH2F("hmm2ppippimVw","Missing Mass2 of p,#pi^{+},#pi^{-} vs W",150, 0, 4, 600,-0.003,0.003));
	ret->Add(new TH2F("hmmppippimVw", "Missing Mass of p,#pi^{+},#pi^{-} vs W", 150, 0, 4, 20000, -0.1, 0.1));
	ret->Add(new TH2F("hmm2ppipVw",   "Missing Mass2 of p,#pi^{+} vs W",        150, 0, 4, 100, -0.02, 0.16));
	ret->Add(new TH2F("hmmppipVw",    "Missing Mass of p,#pi^{+} vs W",         150, 0, 4, 100, 0.00, 0.40));
	ret->Add(new TH2F("hmm2ppimVw",   "Missing Mass2 of p,#pi^{-} vs W",        150, 0, 4, 100, -0.02, 0.16));
	ret->Add(new TH2F("hmmppimVw",    "Missing Mass of p,#pi^{-} vs W",         150, 0, 4, 100, 0.00, 0.40));
	ret->Add(new TH2F("hmm2pippimVw", "Missing Mass2 of #pi^{+},#pi^{-} vs W",  100, 0, 4, 100,  0.8, 1.2));
	ret->Add(new TH2F("hmmpippimVw",  "Missing Mass of #pi^{+},#pi^{-} vs W",   100, 0, 4, 100,  0.8, 1.2));
	return ret;
}

TObjArray* DataAna::makeYields()
{
	Int_t numHists = 3; 
	TObjArray *ret = new TObjArray(numHists);
	Int_t hdim = 8;
	
	bngQ2.bins = 100;
	bngQ2.xmin = 0.0;
	bngQ2.xmax = 5.0;
	
	bngW.bins  = 400;
	bngW.xmin = 1.0;
	bngW.xmax = 3.0;
	
	bngMppip.bins = 80; //~22 MeV/bin
	bngMppip.xmin = 0.938 + 0.140;
	bngMppip.xmax = bngW.xmax - 0.140;
	
	bngMppim.bins = 80; //~22 MeV/bin;
	bngMppim.xmin = 0.938 + 0.140;
	bngMppim.xmax = bngW.xmax - 0.140;
	
	bngMpippim.bins = 80; //~22 MeV/bin;
	bngMpippim.xmin = 0.140 + 0.140;
	bngMpippim.xmax = bngW.xmax - 0.938;
	
	bngTheta.bins = 10; //10;
	bngTheta.xmin = 0;
	bngTheta.xmax = 180;
	
	bngPhi.bins = 10; //10;
	bngPhi.xmin = 0;
	bngPhi.xmax = 360;
	
	bngAlpha.bins = 1;
	bngAlpha.xmin = 0;
	bngAlpha.xmax = 360;
	
	/* Varset 1*/
	//                    {  h, Q2,         W,         Mppip,         Mpippim,         theta_pim,     phi_pim,     alpha[p'pip][ppim]}
	Int_t bins1[]    =    {  3, bngQ2.bins, bngW.bins, bngMppip.bins, bngMpippim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
	Double_t xmin1[] =    { -1, bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMpippim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
	Double_t xmax1[] =    {  2, bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMpippim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
	THnSparse* hN1 = new THnSparseF("yield_varset1", 
	"h, Q^{2}, W, M_{p#pi^{+}}, M_{#pi^{+}#pi^{-}}, #theta_{#pi^{-}}, #phi_{#pi^{-}}, #alpha_{[p^{'}#pi^{+}][p#pi^{-}]}", 
	hdim, bins1, xmin1, xmax1);
	hN1->Sumw2();
	gDirectory->Append(hN1);
	ret->Add(hN1);

	//atrivedi:080213	
	/* Varset 2*/
	//                    {  h, Q2,         W,         Mppip,         Mpippim,         theta_p,       phi_p,       alpha[pippim][p'p]}
	/*Int_t bins2[]    =    {  3, bngQ2.bins, bngW.bins, bngMppip.bins, bngMpippim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
	Double_t xmin2[] =    { -1, bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMpippim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
	Double_t xmax2[] =    {  2, bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMpippim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
	THnSparse* hN2 = new THnSparseF("yield_varset2", 
	"h, Q^{2}, W, M_{p#pi^{+}}, M_{#pi^{+}#pi^{-}}, #theta_{p}, #phi_{p}, #alpha_{[#pi^{+}#pi^{-}][p^{'}p]}", 
	hdim, bins2, xmin2, xmax2);
	hN2->Sumw2();
	gDirectory->Append(hN2);
	ret->Add(hN2);*/
	
	/* Varset 3*/
	//                    {  h,  Q2,         W,         Mppip,         Mppim,         theta_pip,     phi_pip,     alpha[p'pim][ppip]}
	/*Int_t bins3[]    =    {  3,  bngQ2.bins, bngW.bins, bngMppip.bins, bngMppim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
	Double_t xmin3[] =    { -1,  bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMppim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
	Double_t xmax3[] =    {  2,  bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMppim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
	THnSparse* hN3 = new THnSparseF("yield_varset3", 
	"h, Q^{2}, W, M_{p#pi^{+}}, M_{p#pi^{-}}, #theta_{#pi^{+}}, #phi_{#pi^{+}}, #alpha_{[p^{'}#pi^{-}][p#pi^{+}]}", 
	hdim, bins3, xmin3, xmax3);
	hN3->Sumw2();
	gDirectory->Append(hN3);
	ret->Add(hN3);*/
	
	return ret;
}

void DataAna::writeHists(TObjArray** hists, TDirectory *dirout)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if (iSector>0) continue;
		Info("DataAna::writeHists()", "here1");
		Info("DataAna::writeHists()", "cd(sector%d)", iSector);
		dirout->cd(TString::Format("sector%d", iSector));
		Info("DataAna::writeHists()", "here2");
		hists[iSector]->Write();
		Info("DataAna::writeHists()", "here3");
	}
	return;
}

void DataAna::deleteHists(TObjArray** hists)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		Info("DataAna::deleteHists()", "here1");
		//delete hists[iSector];
		Info("DataAna::deleteHists()", "here2");
	}
	return;
}

void DataAna::fillHistsEid(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	//Int_t iSector = dAna->eid.sector;
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (iSector == 0) || (iSector == eid.sector) ) {
						
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			h1->Fill(eid.ec_ei, eid.ec_eo);
					
			TH2* h2 = (TH2*)hists[iSector]->At(1);
			h2->Fill(eid.ec_ei/eid.p, eid.ec_eo/eid.p);
	
			TH2* h3 = (TH2*)hists[iSector]->At(2);
			h3->Fill(eid.p, eid.etot/eid.p);
	
			TH1* h4 = (TH1*)hists[iSector]->At(3);
			h4->Fill(eid.nphe);
	
			TH2* h5 = (TH2*)hists[iSector]->At(4);
			h5->Fill(eid.cc_segm, eid.cc_theta);
	
			TH2* h6 = (TH2*)hists[iSector]->At(5);
			h6->Fill(eid.p, eid.beta);
	
			TH2* h7 = (TH2*)hists[iSector]->At(6);
			h7->Fill(eid.p, eid.dtE);
		}
	}
	/*TH2* h1 = (TH2*)hists->At(0);
	h1->Fill(eid.ec_ei, eid.ec_eo);
	
	TH2* h2 = (TH2*)hists->At(1);
	h2->Fill(eid.ec_ei/eid.p, eid.ec_eo/eid.p);
	
	TH2* h3 = (TH2*)hists->At(2);
	h3->Fill(eid.p, eid.etot/eid.p);
	
	TH1* h4 = (TH1*)hists->At(3);
	h4->Fill(eid.nphe);
	
	TH2* h5 = (TH2*)hists->At(4);
	h5->Fill(eid.cc_segm, eid.cc_theta);
	
	TH2* h6 = (TH2*)hists->At(5);
	h6->Fill(eid.p, eid.beta);
	
	TH2* h7 = (TH2*)hists->At(6);
	h7->Fill(eid.p, eid.dtE);*/
}

void DataAna::fillHistsEFid(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (iSector == 0) || (iSector == efid.sector) ) {
		
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			h1->Fill(efid.theta, efid.phi);
			
			TH2* h2 = (TH2*)hists[iSector]->At(1);
			h2->Fill(efid.dc_xsc, efid.dc_ysc);
			
			TH2* h3 = (TH2*)hists[iSector]->At(2);
			h3->Fill(efid.ech_x, efid.ech_y);
		}
	}
}

void DataAna::fillHistsMomCor(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (iSector == 0) || (iSector == mom.sector) ) {
		
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			h1->Fill(mom.p, mom.dp);
	
			TH1* h2 = (TH1*)hists[iSector]->At(1);
			h2->Fill(mom.dcx);
	
			TH1* h3 = (TH1*)hists[iSector]->At(2);
			h3->Fill(mom.dcy);
	
			TH1* h4 = (TH3*)hists[iSector]->At(3);
			h4->Fill(mom.dcz);
	
			TH1* h5 = (TH1*)hists[iSector]->At(4);
			h5->Fill(mom.dp);
		}
	}
}

void DataAna::fillHistsPid(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	Int_t hIdx = 0;
	if(pid.h10IdxP > 0){
		hIdx = 0;
		for(Int_t iSector=0;iSector<NSECTORS;iSector++){
				if ( (iSector == 0) || (iSector == pid.sectorP) ) {
			
				TH2* h1 = (TH2*)hists[iSector]->At(hIdx);
				h1->Fill(pid.pP, pid.betaP);
				TH2* h2 = (TH2*)hists[iSector]->At(hIdx+1);
				h2->Fill(pid.pP, pid.dtP);
		
				TH2* h3 = (TH2*)hists[iSector]->At(hIdx+2);
				h3->Fill(pid.P_ec_ei, pid.P_ec_eo);
				TH2* h4 = (TH2*)hists[iSector]->At(hIdx+3);
				h4->Fill(pid.P_ec_ei/pid.pP, pid.P_ec_eo/pid.pP);
				TH2* h5 = (TH2*)hists[iSector]->At(hIdx+4);
				h5->Fill(pid.pP, pid.P_etot/pid.pP);
				TH1* h6 = (TH1*)hists[iSector]->At(hIdx+5);
				h6->Fill(pid.P_nphe);
				TH2* h7 = (TH2*)hists[iSector]->At(hIdx+6);
				h7->Fill(pid.P_cc_segm, pid.P_cc_theta);
			}
		}
	}
	if(pid.h10IdxPip > 0){
		hIdx = 7;
		for(Int_t iSector=0;iSector<NSECTORS;iSector++){
			if ( (iSector == 0) || (iSector == pid.sectorPip) ) {
			
				TH2* h1 = (TH2*)hists[iSector]->At(hIdx);
				h1->Fill(pid.pPip, pid.betaPip);
				TH2* h2 = (TH2*)hists[iSector]->At(hIdx+1);
				h2->Fill(pid.pPip, pid.dtPip);
		
				TH2* h3 = (TH2*)hists[iSector]->At(hIdx+2);
				h3->Fill(pid.Pip_ec_ei, pid.Pip_ec_eo);
				TH2* h4 = (TH2*)hists[iSector]->At(hIdx+3);
				h4->Fill(pid.Pip_ec_ei/pid.pPip, pid.Pip_ec_eo/pid.pPip);
				TH2* h5 = (TH2*)hists[iSector]->At(hIdx+4);
				h5->Fill(pid.pPip, pid.Pip_etot/pid.pPip);
				TH1* h6 = (TH1*)hists[iSector]->At(hIdx+5);
				h6->Fill(pid.Pip_nphe);
				TH2* h7 = (TH2*)hists[iSector]->At(hIdx+6);
				h7->Fill(pid.Pip_cc_segm, pid.Pip_cc_theta);
			}
		}
	}
	if(pid.h10IdxPim > 0){
		hIdx = 14;
		for(Int_t iSector=0;iSector<NSECTORS;iSector++){
			if ( (iSector == 0) || (iSector == pid.sectorPim) ) {
			
				TH2* h1 = (TH2*)hists[iSector]->At(hIdx);
				h1->Fill(pid.pPim, pid.betaPim);
				TH2* h2 = (TH2*)hists[iSector]->At(hIdx+1);
				h2->Fill(pid.pPim, pid.dtPim);
		
				TH2* h3 = (TH2*)hists[iSector]->At(hIdx+2);
				h3->Fill(pid.Pim_ec_ei, pid.Pim_ec_eo);
				TH2* h4 = (TH2*)hists[iSector]->At(hIdx+3);
				h4->Fill(pid.Pim_ec_ei/pid.pPim, pid.Pim_ec_eo/pid.pPim);
				TH2* h5 = (TH2*)hists[iSector]->At(hIdx+4);
				h5->Fill(pid.pPim, pid.Pim_etot/pid.pPim);
				TH1* h6 = (TH1*)hists[iSector]->At(hIdx+5);
				h6->Fill(pid.Pim_nphe);
				TH2* h7 = (TH2*)hists[iSector]->At(hIdx+6);
				h7->Fill(pid.Pim_cc_segm, pid.Pim_cc_theta);
			}
		}
	}
}

void DataAna::fillHistsEkin(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	DataEkin *ekin = &eKin;
	if (useMc) ekin = &eKin_mc;
	
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (iSector != 0) && (iSector != ekin->sector) ) continue;
		
		TH1* h1 = (TH1*)hists[iSector]->At(0);
		h1->Fill(ekin->W);
	
		TH2* h2 = (TH2*)hists[iSector]->At(1);
		h2->Fill(ekin->W,ekin->Q2);
	
		TH1* h3 = (TH1*)hists[iSector]->At(2);
		h3->Fill(ekin->phi);
	
		TH2* h4 = (TH2*)hists[iSector]->At(3);
		h4->Fill(ekin->theta1,ekin->E1);
	
		TH2* h5 = (TH2*)hists[iSector]->At(4);
		h5->Fill(ekin->theta1,ekin->Q2);
	
		TH2* h6 = (TH2*)hists[iSector]->At(5);
		h6->Fill(ekin->E1, ekin->Q2);
	}
}

void DataAna::fillHistsEkin(TObjArray* hists, Bool_t useMc /* = kFALSE */)
{
	DataEkin *ekin = &eKin;
	if (useMc) ekin = &eKin_mc;
	
	TH1* h1 = (TH1*)hists->At(0);
	h1->Fill(ekin->W);
	
	TH2* h2 = (TH2*)hists->At(1);
	h2->Fill(ekin->W,ekin->Q2);
	
	TH1* h3 = (TH1*)hists->At(2);
	h3->Fill(ekin->phi);
	
	TH2* h4 = (TH2*)hists->At(3);
	h4->Fill(ekin->theta1,ekin->E1);
	
	TH2* h5 = (TH2*)hists->At(4);
	h5->Fill(ekin->theta1,ekin->Q2);
	
	TH2* h6 = (TH2*)hists->At(5);
	h6->Fill(ekin->E1, ekin->Q2);
	
}

void DataAna::fillHistsMM(TObjArray *hists, Bool_t useMc /* = kFALSE */)
{
	Int_t hIdx = 0;
	Data2pi *tp = &d2pi;
	if (useMc) tp = &d2pi_mc;
		
	if ((useMc || (h10idxP>0 && h10idxPip>0 && h10idxPim>0)) ) {	
		hIdx = 0;
		TH2* h1 = (TH2*)hists->At(hIdx);
		h1->Fill(tp->W,tp->mm2ppippim);
		
		TH2* h2 = (TH2*)hists->At(hIdx+1);
		h2->Fill(tp->W,tp->mmppippim);
	}
	
    if ((useMc || (h10idxP>0 && h10idxPip>0 && h10idxPim==-1)) ) {
		hIdx = 2;
		TH2* h1 = (TH2*)hists->At(hIdx);
		h1->Fill(tp->W,tp->mm2ppip);
		
		TH2* h2 = (TH2*)hists->At(hIdx+1);
		h2->Fill(tp->W,tp->mmppip);
	}
	
	if ((useMc || (h10idxP>0 && h10idxPim>0 && h10idxPip==-1)) ) {
		hIdx = 4;
		TH2* h1 = (TH2*)hists->At(hIdx);
		h1->Fill(tp->W,tp->mm2ppim);
		
		TH2* h2 = (TH2*)hists->At(hIdx+1);
		h2->Fill(tp->W,tp->mmppim);
	}
	
	if ((useMc || (h10idxPip>0 && h10idxPim>0 && h10idxP==-1)) ) {
		hIdx = 6;
		TH2* h1 = (TH2*)hists->At(hIdx);
		h1->Fill(tp->W,tp->mm2pippim);
		
		TH2* h2 = (TH2*)hists->At(hIdx+1);
		h2->Fill(tp->W,tp->mmpippim);
	}
	
}

void DataAna::fillYields(TObjArray *hists, Bool_t useMc /* = kFALSE */)
{
	Data2pi *tp = &d2pi;
	if (useMc) tp = &d2pi_mc;
	THnSparse* hN1 = (THnSparse*)hists->At(0);
	Double_t coord1[] = { tp->h, tp->Q2, tp->W, tp->varset1.M1, tp->varset1.M2, tp->varset1.theta, tp->varset1.phi, tp->varset1.alpha };
	hN1->Fill(coord1);
	
	//atrivedi:080231
	/*THnSparse* hN2 = (THnSparse*)hists->At(1);
	Double_t coord2[] = { tp->h, tp->Q2, tp->W, tp->varset2.M1, tp->varset2.M2, tp->varset2.theta, tp->varset2.phi, tp->varset2.alpha  };
	hN2->Fill(coord2);
	
	THnSparse* hN3 = (THnSparse*)hists->At(2);
	Double_t coord3[] = { tp->h, tp->Q2, tp->W, tp->varset3.M1, tp->varset3.M2, tp->varset3.theta, tp->varset3.phi, tp->varset3.alpha  };
	hN3->Fill(coord3);*/
		
}
