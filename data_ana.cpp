#include "data_ana.h"
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TMath.h>
#include "particle_constants.h"
//#include "q2w_bng.h"
#include "h8_bng.h"

#include <iostream>
using namespace std;

using namespace TMath;
using namespace AnalysisConstants;
using namespace ParticleConstants;

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
	pfid.Clear();
	pfid_elast.Clear();
	eeff.Clear();
	skimq.Clear();
	skimq_elast.Clear();
	mom.Clear();
	pid.Clear();
	pidnew.Clear();
	peff.Clear();
	pid_elast.Clear();
	pid_elastnew.Clear();
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

void DataAna::makeHistsPFid(TObjArray** hists, TDirectory* dirout)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(3);
		
		if (iSector==0)	{
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 0)",     100,0,60, 100,0,360));
			hists[iSector]->Add(new TH2F("h_pip_phiVtheta","#phi vs. #theta for pip(sector 0)", 100,0,120,100,0,360));
			hists[iSector]->Add(new TH2F("h_pim_phiVtheta","#phi vs. #theta for pim(sector 0)", 100,0,120,100,0,360));
		}
		if (iSector==1) {
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 1)",     100,0,60, 100,0,360));
			hists[iSector]->Add(new TH2F("h_pip_phiVtheta","#phi vs. #theta for pip(sector 1)", 100,0,120,100,0,360));
			hists[iSector]->Add(new TH2F("h_pim_phiVtheta","#phi vs. #theta for pim(sector 1)", 100,0,120,100,0,360));
		}
		if (iSector==2) {
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 2)",     100,0,60, 100,0,360));
			hists[iSector]->Add(new TH2F("h_pip_phiVtheta","#phi vs. #theta for pip(sector 2)", 100,0,120,100,0,360));
			hists[iSector]->Add(new TH2F("h_pim_phiVtheta","#phi vs. #theta for pim(sector 2)", 100,0,120,100,0,360));
		}
		if (iSector==3)	{
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 3)",     100,0,60, 100,0,360));
			hists[iSector]->Add(new TH2F("h_pip_phiVtheta","#phi vs. #theta for pip(sector 3)", 100,0,120,100,0,360));
			hists[iSector]->Add(new TH2F("h_pim_phiVtheta","#phi vs. #theta for pim(sector 3)", 100,0,120,100,0,360));
		}
		if (iSector==4) {
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 4)",     100,0,60, 100,0,360));
			hists[iSector]->Add(new TH2F("h_pip_phiVtheta","#phi vs. #theta for pip(sector 4)", 100,0,120,100,0,360));
			hists[iSector]->Add(new TH2F("h_pim_phiVtheta","#phi vs. #theta for pim(sector 4)", 100,0,120,100,0,360));
		}
		if (iSector==5) {
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 5)",     100,0,60, 100,0,360));
			hists[iSector]->Add(new TH2F("h_pip_phiVtheta","#phi vs. #theta for pip(sector 5)", 100,0,120,100,0,360));
			hists[iSector]->Add(new TH2F("h_pim_phiVtheta","#phi vs. #theta for pim(sector 5)", 100,0,120,100,0,360));
		}
		if (iSector==6) {
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 6)",     100,0,60, 100,0,360));
			hists[iSector]->Add(new TH2F("h_pip_phiVtheta","#phi vs. #theta for pip(sector 6)", 100,0,120,100,0,360));
			hists[iSector]->Add(new TH2F("h_pim_phiVtheta","#phi vs. #theta for pim(sector 6)", 100,0,120,100,0,360));
		}
	}
	//return ret;
	return;
}

void DataAna::makeHistsPFidElast(TObjArray** hists, TDirectory* dirout)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(1);
		
		if (iSector==0)	{
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 0)",     100,0,60, 100,0,360));
		}
		if (iSector==1) {
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 1)",     100,0,60, 100,0,360));
		}
		if (iSector==2) {
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 2)",     100,0,60, 100,0,360));
		}
		if (iSector==3)	{
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 3)",     100,0,60, 100,0,360));
		}
		if (iSector==4) {
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 4)",     100,0,60, 100,0,360));
		}
		if (iSector==5) {
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 5)",     100,0,60, 100,0,360));
		}
		if (iSector==6) {
			hists[iSector]->Add(new TH2F("h_p_phiVtheta","#phi vs. #theta for p(sector 6)",     100,0,60, 100,0,360));
		}
	}
	//return ret;
	return;
}

void DataAna::makeHistsEEff(TObjArray** hists, TDirectory* dirout)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(1);
		
		if (iSector==0)	hists[iSector]->Add(new TH2F("h_e_thetaVp","e #theta vs. p(sector 0)", 100,1,5,100,0,60));
		if (iSector==1) hists[iSector]->Add(new TH2F("h_e_thetaVp","e #theta vs. p(sector 1)", 100,1,5,100,0,60));
		if (iSector==2) hists[iSector]->Add(new TH2F("h_e_thetaVp","e #theta vs. p(sector 2)", 100,1,5,100,0,60));
		if (iSector==3)	hists[iSector]->Add(new TH2F("h_e_thetaVp","e #theta vs. p(sector 3)", 100,1,5,100,0,60));
		if (iSector==4) hists[iSector]->Add(new TH2F("h_e_thetaVp","e #theta vs. p(sector 4)", 100,1,5,100,0,60));
		if (iSector==5) hists[iSector]->Add(new TH2F("h_e_thetaVp","e #theta vs. p(sector 5)", 100,1,5,100,0,60));
		if (iSector==6) hists[iSector]->Add(new TH2F("h_e_thetaVp","e #theta vs. p(sector 6)", 100,1,5,100,0,60));
	}
	//return ret;
	return;
}

void DataAna::makeHistsPEff(TObjArray** hists, TDirectory* dirout)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(3);
		
		if (iSector==0)	{
			hists[iSector]->Add(new TH2F("h_p_thetaVp",  "#theta vs. p for p(sector 0)",      100,0,4,100,0,60));
			hists[iSector]->Add(new TH2F("h_pip_thetaVp","#theta vs. p for #pi^{+}(sector 0)",100,0,3,100,0,120));
			hists[iSector]->Add(new TH2F("h_pim_thetaVp","#theta vs. p for #pi^{-}(sector 0)",100,0,3,100,0,120));
		}
		if (iSector==1) {
			hists[iSector]->Add(new TH2F("h_p_thetaVp",  "#theta vs. p for p(sector 1)",      100,0,4,100,0,60));
			hists[iSector]->Add(new TH2F("h_pip_thetaVp","#theta vs. p for #pi^{+}(sector 1)",100,0,3,100,0,120));
			hists[iSector]->Add(new TH2F("h_pim_thetaVp","#theta vs. p for #pi^{-}(sector 1)",100,0,3,100,0,120));
		}
		if (iSector==2) {
			hists[iSector]->Add(new TH2F("h_p_thetaVp",  "#theta vs. p for p(sector 2)",      100,0,4,100,0,60));
			hists[iSector]->Add(new TH2F("h_pip_thetaVp","#theta vs. p for #pi^{+}(sector 2)",100,0,3,100,0,120));
			hists[iSector]->Add(new TH2F("h_pim_thetaVp","#theta vs. p for #pi^{-}(sector 2)",100,0,3,100,0,120));
		}
		if (iSector==3)	{
			hists[iSector]->Add(new TH2F("h_p_thetaVp",  "#theta vs. p for p(sector 3)",      100,0,4,100,0,60));
			hists[iSector]->Add(new TH2F("h_pip_thetaVp","#theta vs. p for #pi^{+}(sector 3)",100,0,3,100,0,120));
			hists[iSector]->Add(new TH2F("h_pim_thetaVp","#theta vs. p for #pi^{-}(sector 3)",100,0,3,100,0,120));
		}
		if (iSector==4) {
			hists[iSector]->Add(new TH2F("h_p_thetaVp",  "#theta vs. p for p(sector 4)",      100,0,4,100,0,60));
			hists[iSector]->Add(new TH2F("h_pip_thetaVp","#theta vs. p for #pi^{+}(sector 4)",100,0,3,100,0,120));
			hists[iSector]->Add(new TH2F("h_pim_thetaVp","#theta vs. p for #pi^{-}(sector 4)",100,0,3,100,0,120));
		}
		if (iSector==5) {
			hists[iSector]->Add(new TH2F("h_p_thetaVp",  "#theta vs. p for p(sector 5)",      100,0,4,100,0,60));
			hists[iSector]->Add(new TH2F("h_pip_thetaVp","#theta vs. p for #pi^{+}(sector 5)",100,0,3,100,0,120));
			hists[iSector]->Add(new TH2F("h_pim_thetaVp","#theta vs. p for #pi^{-}(sector 5)",100,0,3,100,0,120));
		}
		if (iSector==6) {
			hists[iSector]->Add(new TH2F("h_p_thetaVp",  "#theta vs. p for p(sector 6)",      100,0,4,100,0,60));
			hists[iSector]->Add(new TH2F("h_pip_thetaVp","#theta vs. p for #pi^{+}(sector 6)",100,0,3,100,0,120));
			hists[iSector]->Add(new TH2F("h_pim_thetaVp","#theta vs. p for #pi^{-}(sector 6)",100,0,3,100,0,120));
		}
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
		hists[iSector]->Add(new TH2F("hP_dtVp",TString::Format("#Deltat vs p for idtfd proton(sector%d)", iSector),100,0,5,100,-2,2));
		hists[iSector]->Add(new TH2F("hP_eoutVein",TString::Format("#DeltaE_{out} vs. #DeltaE_{in} for idtfd p(sector%d)", iSector),150, 0, 1, 150, 0, 1));
		hists[iSector]->Add(new TH2F("hP_sfoutVsfin",TString::Format("SF_{out} vs. SF_{in} for idtfd p(sector%d)", iSector),100, 0, 0.5, 100, 0, 0.5));
		hists[iSector]->Add(new TH2F("hP_sftotVp",TString::Format("SF_{tot} vs. p for idtfd p(sector%d)", iSector), 160, 0, 5, 100, 0, 0.5));
		hists[iSector]->Add(new TH1F("hP_nphe", TString::Format("nphe for idtfd p(sector%d)", iSector), 100, 0, 300));
		hists[iSector]->Add(new TH2F("hP_CCthetaVCCseg",TString::Format("CC_{#theta} vs. CC_{seg} for idtf p(sector%d)", iSector), 20, 0, 20, 100, 0, 50));
		
		hists[iSector]->Add(new TH2F("hPip_betaVp", TString::Format("#beta vs p for idtfd #pi^{+}(sector%d)", iSector),100, 0, 4, 100, 0, 1.5));
		hists[iSector]->Add(new TH2F("hPip_dtVp", TString::Format("#Deltat vs p for idtfd #pi^{+}(sector%d)", iSector),100,0,5,100,-2,2));
		hists[iSector]->Add(new TH2F("hPip_eoutVein", TString::Format("#DeltaE_{out} vs. #DeltaE_{in} for idtfd #pi^{+}(sector%d)", iSector),150, 0, 1, 150, 0, 1));
		hists[iSector]->Add(new TH2F("hPip_sfoutVsfin", TString::Format("SF_{out} vs. SF_{in} for idtfd #pi^{+}(sector%d)", iSector),100, 0, 0.5, 100, 0, 0.5));
		hists[iSector]->Add(new TH2F("hPip_sftotVp", TString::Format("SF_{tot} vs. p for idtfd #pi^{+}(sector%d)", iSector), 160, 0, 5, 100, 0, 0.5));
		hists[iSector]->Add(new TH1F("hPip_nphe", TString::Format("nphe for idtfd #pi^{+}(sector%d)", iSector), 100, 0, 300));
		hists[iSector]->Add(new TH2F("hPip_CCthetaVCCseg", TString::Format("CC_{#theta} vs. CC_{seg} for idtf #pi^{+}(sector%d)", iSector), 20, 0, 20, 100, 0, 50));
		
		hists[iSector]->Add(new TH2F("hPim_betaVp", TString::Format("#beta vs p for idtfd #pi^{-}(sector%d)", iSector),100, 0, 4, 100, 0, 1.5));
		hists[iSector]->Add(new TH2F("hPim_dtVp", TString::Format("#Deltat vs p for idtfd #pi^{-}(sector%d)", iSector),100,0,5,100,-2,2));
		hists[iSector]->Add(new TH2F("hPim_eoutVein", TString::Format("#DeltaE_{out} vs. #DeltaE_{in} for idtfd #pi^{-}(sector%d)", iSector),150, 0, 1, 150, 0, 1));
		hists[iSector]->Add(new TH2F("hPim_sfoutVsfin", TString::Format("SF_{out} vs. SF_{in} for idtfd #pi^{-}(sector%d)", iSector),100, 0, 0.5, 100, 0, 0.5));
		hists[iSector]->Add(new TH2F("hPim_sftotVp", TString::Format("SF_{tot} vs. p for idtfd #pi^{-}(sector%d)", iSector), 160, 0, 5, 100, 0, 0.5));
		hists[iSector]->Add(new TH1F("hPim_nphe", TString::Format("nphe for idtfd #pi^{-}(sector%d)", iSector), 100, 0, 300));
		hists[iSector]->Add(new TH2F("hPim_CCthetaVCCseg", TString::Format("CC_{#theta} vs. CC_{seg} for idtf #pi^{-}(sector%d)", iSector), 20, 0, 20, 100, 0, 50));
	}
	//return ret;
	return;
}

void DataAna::makeHistsPidMon(TObjArray** hists, TDirectory* dirout)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(5);
	
		hists[iSector]->Add(new TH2F("h_betaVp_pos",TString::Format("#beta v. p for +ve particles (sector%d)", iSector),100,0,5,100,0,1.2));
		hists[iSector]->Add(new TH2F("h_dtVp_pos_p",TString::Format("#Deltat v. p for +ve particles under p hyp.(sector%d)", iSector),100,0,5,250,-5,5));
		hists[iSector]->Add(new TH2F("h_dtVp_pos_pip",TString::Format("#Deltat v. p for +ve particles under #pi^{+} hyp.(sector%d)", iSector),100,0,5,250,-5,5));
				
		hists[iSector]->Add(new TH2F("h_betaVp_neg",TString::Format("#beta v. p for -ve(sector%d)", iSector),100,0,5,100,0,1.2));
		hists[iSector]->Add(new TH2F("h_dtVp_neg_pim",TString::Format("#Deltat v. p for -ve particles under #pi^{-} hyp.(sector%d)", iSector),100,0,5,250,-5,5));
		
	}
	return;
}

void DataAna::makeHistsPidCut(TObjArray** hists, TDirectory* dirout)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(6);
	
		hists[iSector]->Add(new TH2F("h_betaVp_p",TString::Format("#beta v. p for +ve particles idtfd as p (sector%d)", iSector),100,0,5,100,0,1.2));
		hists[iSector]->Add(new TH2F("h_dtVp_p",TString::Format("#Deltat v. p for +ve particles idtfd as p (sector%d)", iSector),100,0,5,250,-5,5));

		hists[iSector]->Add(new TH2F("h_betaVp_pip",TString::Format("#beta v. p for +ve particles idtfd as #pi^{+} (sector%d)", iSector),100,0,5,100,0,1.2));
		hists[iSector]->Add(new TH2F("h_dtVp_pip",TString::Format("#Deltat v. p for +ve particles idtfd as #pi^{+} (sector%d)", iSector),100,0,5,250,-5,5));
				
		hists[iSector]->Add(new TH2F("h_betaVp_pim",TString::Format("#beta v. p for -ve particles idtfd as #pi^{-}(sector%d)", iSector),100,0,5,100,0,1.2));
		hists[iSector]->Add(new TH2F("h_dtVp_pim",TString::Format("#Deltat v. p for -ve particles idtfd as #pi^{-} (sector%d)", iSector),100,0,5,250,-5,5));
		
	}
	return;
}

void DataAna::makeHistsPidElast(TObjArray** hists, TDirectory* dirout)
{
	//TObjArray *ret = new TObjArray(21);
	TObjArray* ret[7];
	//TDirectory* dirpid = dirout->mkdir("pid");
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(7);
	
		hists[iSector]->Add(new TH2F("hP_betaVp",TString::Format("#beta vs p for idtfd proton(sector%d)", iSector),100, 0, 4, 100, 0, 1.5));
		hists[iSector]->Add(new TH2F("hP_dtVp",TString::Format("#Deltat vs p for idtfd proton(sector%d)", iSector),100,0,5,100,-2,2));
		hists[iSector]->Add(new TH2F("hP_eoutVein",TString::Format("#DeltaE_{out} vs. #DeltaE_{in} for idtfd p(sector%d)", iSector),150, 0, 1, 150, 0, 1));
		hists[iSector]->Add(new TH2F("hP_sfoutVsfin",TString::Format("SF_{out} vs. SF_{in} for idtfd p(sector%d)", iSector),100, 0, 0.5, 100, 0, 0.5));
		hists[iSector]->Add(new TH2F("hP_sftotVp",TString::Format("SF_{tot} vs. p for idtfd p(sector%d)", iSector), 160, 0, 5, 100, 0, 0.5));
		hists[iSector]->Add(new TH1F("hP_nphe", TString::Format("nphe for idtfd p(sector%d)", iSector), 100, 0, 300));
		hists[iSector]->Add(new TH2F("hP_CCthetaVCCseg",TString::Format("CC_{#theta} vs. CC_{seg} for idtf p(sector%d)", iSector), 20, 0, 20, 100, 0, 50));
	}
	//return ret;
	return;
}

void DataAna::makeHistsPidElastMon(TObjArray** hists, TDirectory* dirout)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(3);
	
		hists[iSector]->Add(new TH2F("h_betaVp_pos",TString::Format("#beta v. p for +ve particles (sector%d)", iSector),100,0,5,100,0,1.2));
		hists[iSector]->Add(new TH2F("h_dtVp_pos_p",TString::Format("#Deltat v. p for +ve particles under p hyp.(sector%d)", iSector),100,0,5,250,-5,5));
		hists[iSector]->Add(new TH2F("h_dtVp_pos_pip",TString::Format("#Deltat v. p for +ve particles under #pi^{+} hyp.(sector%d)", iSector),100,0,5,250,-5,5));
	}
	return;
}

void DataAna::makeHistsPidElastCut(TObjArray** hists, TDirectory* dirout)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if ( (dirout->GetDirectory(TString::Format("sector%d",iSector))) == NULL )
			dirout->mkdir(TString::Format("sector%d",iSector))->cd();
			
		hists[iSector] = new TObjArray(2);
	
		hists[iSector]->Add(new TH2F("h_betaVp_p",TString::Format("#beta v. p for +ve particles idtfd as p (sector%d)", iSector),100,0,5,100,0,1.2));
		hists[iSector]->Add(new TH2F("h_dtVp_p",TString::Format("#Deltat v. p for +ve particles idtfd as p (sector%d)", iSector),100,0,5,250,-5,5));
	}
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
	TObjArray *ret = new TObjArray(16);
	ret->Add(new TH1F("hmm2ppippim","Missing Mass2 of p,#pi^{+},#pi^{-}",600,   -0.003,0.003));
	ret->Add(new TH1F("hmmppippim", "Missing Mass of p,#pi^{+},#pi^{-}", 20000, -0.10,  0.1));
	ret->Add(new TH1F("hmm2ppip",   "Missing Mass2 of p,#pi^{+}",        100,   -0.02, 1.00));//0.16));
	ret->Add(new TH1F("hmmppip",    "Missing Mass of p,#pi^{+}",         100,    0.00, 1.00));//0.40));
	ret->Add(new TH1F("hmm2ppim",   "Missing Mass2 of p,#pi^{-}",        100,   -0.02, 1.00));//0.16));
	ret->Add(new TH1F("hmmppim",    "Missing Mass of p,#pi^{-}",         100,    0.00, 1.00));//0.40));
	ret->Add(new TH1F("hmm2pippim", "Missing Mass2 of #pi^{+},#pi^{-}",  100,    0.80,  1.2));
	ret->Add(new TH1F("hmmpippim",  "Missing Mass of #pi^{+},#pi^{-}",   100,    0.80,  1.2));
	ret->Add(new TH2F("hmm2ppippimVw","Missing Mass2 of p,#pi^{+},#pi^{-} vs W",150, 0, 3, 600,   -0.003,0.003));
	ret->Add(new TH2F("hmmppippimVw", "Missing Mass of p,#pi^{+},#pi^{-} vs W", 150, 0, 3, 20000, -0.10,  0.1));
	ret->Add(new TH2F("hmm2ppipVw",   "Missing Mass2 of p,#pi^{+} vs W",        150, 0, 3, 100,   -0.02, 1.00));//0.16));
	ret->Add(new TH2F("hmmppipVw",    "Missing Mass of p,#pi^{+} vs W",         150, 0, 3, 100,    0.00, 1.00));//0.40));
	ret->Add(new TH2F("hmm2ppimVw",   "Missing Mass2 of p,#pi^{-} vs W",        150, 0, 3, 100,   -0.02, 1.00));//0.16));
	ret->Add(new TH2F("hmmppimVw",    "Missing Mass of p,#pi^{-} vs W",         150, 0, 3, 100,    0.00, 1.00));//0.40));
	ret->Add(new TH2F("hmm2pippimVw", "Missing Mass2 of #pi^{+},#pi^{-} vs W",  100, 0, 3, 100,    0.80,  1.2));
	ret->Add(new TH2F("hmmpippimVw",  "Missing Mass of #pi^{+},#pi^{-} vs W",   100, 0, 3, 100,    0.80,  1.2));
	return ret;
}

TObjArray* DataAna::makeHistsMMElastic()
{
	TObjArray *ret = new TObjArray(2);
	ret->Add(new TH2F("hmmepVw","Missing Mass of ep vs. W",150,0,3,100,-0.02,1.00));
	ret->Add(new TH1F("hW", "W", 150,0,3));
	return ret;
}

TObjArray** DataAna::makeYields()
{
	TObjArray** ret=new TObjArray* [NBINS_WCRS];
	Int_t hdim=8;
	
	bngQ2.bins=NBINS_Q2; 
	bngQ2.xmin=XMIN_Q2;
	bngQ2.xmax=XMAX_Q2;

	bngTheta.bins=NBINS_THETA;
	bngTheta.xmin=XMIN_THETA;
	bngTheta.xmax=XMAX_THETA;
	
	bngPhi.bins=NBINS_PHI;
	bngPhi.xmin=XMIN_PHI;
	bngPhi.xmax=XMAX_PHI;
	
	bngAlpha.bins=NBINS_ALPHA; 
	bngAlpha.xmin=XMIN_ALPHA;
	bngAlpha.xmax=XMAX_ALPHA;

	for(int i=0;i<NBINS_WCRS;i++){
		TObjArray* h8=new TObjArray(NVARSETS);
	
		bngW.bins=WCRSBIN[i].nbins;
		bngW.xmin=WCRSBIN[i].xmin;
		bngW.xmax=WCRSBIN[i].xmax;
		//bngW.bins=(int)((bngW.xmax-bngW.xmin)/BINW_W);
	
		bngMppip.bins=NBINS_M;
		bngMppip.xmin=MASS_P+MASS_PIP;
		bngMppip.xmax=bngW.xmax-MASS_PIM;
	
		bngMppim.bins=NBINS_M;
		bngMppim.xmin=MASS_P+MASS_PIM;
		bngMppim.xmax=bngW.xmax-MASS_PIM;
	
		bngMpippim.bins=NBINS_M;
		bngMpippim.xmin=MASS_PIP+MASS_PIM;
		bngMpippim.xmax=bngW.xmax-MASS_P;
	
		//! Varsets as per "three assignments(i,ii,iii)"  listed on page 4 of JM06_model_PhysRevC.80.045212.pdf
		/* Varset 1*/
		//                    {  h, Q2,         W,         Mppip,         Mpippim,         theta_pim,     phi_pim,     alpha[p'pip][ppim]}
		Int_t bins1[]    =    {  3, bngQ2.bins, bngW.bins, bngMppip.bins, bngMpippim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
		Double_t xmin1[] =    { -1, bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMpippim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
		Double_t xmax1[] =    {  2, bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMpippim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
		THnSparse* h8_1 = new THnSparseF(TString::Format("h8_%d_%d",i+1,1), 
		"h, Q^{2}, W, M_{p#pi^{+}}, M_{#pi^{+}#pi^{-}}, #theta_{#pi^{-}}, #phi_{#pi^{-}}, #alpha_{[p^{'}#pi^{+}][p#pi^{-}]}", 
		hdim, bins1, xmin1, xmax1);
		h8_1->Sumw2();
		//! Make variable Q2-binning
		//h8_1->GetAxis(1)->Set(kQ2_NAnaBins,kQ2_AnaBins);
		gDirectory->Append(h8_1);
		h8->Add(h8_1);

		/* Varset 2*/
		//                    {  h, Q2,         W,         Mppip,         Mpippim,         theta_p,       phi_p,       alpha[pippim][p'p]}
		Int_t bins2[]    =    {  3, bngQ2.bins, bngW.bins, bngMppip.bins, bngMpippim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
		Double_t xmin2[] =    { -1, bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMpippim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
		Double_t xmax2[] =    {  2, bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMpippim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
		THnSparse* h8_2 = new THnSparseF(TString::Format("h8_%d_%d",i+1,2), 
		"h, Q^{2}, W, M_{p#pi^{+}}, M_{#pi^{+}#pi^{-}}, #theta_{p}, #phi_{p}, #alpha_{[#pi^{+}#pi^{-}][pp^{'}]}", 
		hdim, bins2, xmin2, xmax2);
		h8_2->Sumw2();
		//! Make variable Q2-binning
		//h8_2->GetAxis(1)->Set(kQ2_NAnaBins,kQ2_AnaBins);
		gDirectory->Append(h8_2);
		h8->Add(h8_2);
	
		/* Varset 3*/
		//                    {  h,  Q2,         W,         Mppip,         Mppim,         theta_pip,     phi_pip,     alpha[p'pim][ppip]}
		Int_t bins3[]    =    {  3,  bngQ2.bins, bngW.bins, bngMppip.bins, bngMppim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
		Double_t xmin3[] =    { -1,  bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMppim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
		Double_t xmax3[] =    {  2,  bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMppim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
		THnSparse* h8_3 = new THnSparseF(TString::Format("h8_%d_%d",i+1,3), 
		"h, Q^{2}, W, M_{p#pi^{+}}, M_{p#pi^{-}}, #theta_{#pi^{+}}, #phi_{#pi^{+}}, #alpha_{[p^{'}#pi^{-}][p#pi^{+}]}", 
		hdim, bins3, xmin3, xmax3);
		h8_3->Sumw2();
		//! Make variable Q2-binning
		//h8_3->GetAxis(1)->Set(kQ2_NAnaBins,kQ2_AnaBins);
		gDirectory->Append(h8_3);
		h8->Add(h8_3);

		ret[i]=h8;
	}	
	return ret;
}

TObjArray* DataAna::makeYieldsElastic()
{
	TObjArray* ret=new TObjArray(1);
	Int_t hdim=2;

	Int_t theta_nbins=46;// !1deg/bin (make consistent with Gleb) old: 40; //2 deg/bin
	Int_t theta_min=14;//!efid boundary ~ 13.2 degrees            old: 0;
	Int_t theta_max=60;//!                                        old: 80 

	Int_t phi_nbins=180; //2 deg/bin
	Int_t phi_min=0;
	Int_t phi_max=360;


	//                   {theta,       phi      }
	Int_t bins[]    =    {theta_nbins, phi_nbins};
	Double_t xmin[] =    {theta_min,   phi_min};
	Double_t xmax[] =    {theta_max,   phi_max};
	THnSparse* h2=new THnSparseF("yield","theta, phi",hdim, bins, xmin,xmax);
	h2->Sumw2();
    gDirectory->Append(h2);
	ret->Add(h2);

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

void DataAna::fillHistsPFid(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if (iSector==0) {
			if (pid.h10IdxP>0){
				TH2* h1 = (TH2*)hists[iSector]->At(0);
				float phi=pfid.phi_p+180; //as per fid. cuts implemented by Isupov
				h1->Fill(pfid.theta_p,phi);
			}
			if (pid.h10IdxPip>0){
				TH2* h2 = (TH2*)hists[iSector]->At(1);
				float phi=pfid.phi_pip+180; //as per fid. cuts implemented by Isupov
				h2->Fill(pfid.theta_pip,phi);
				
			}
			if (pid.h10IdxPim>0 ){
				TH2* h3 = (TH2*)hists[iSector]->At(2);
				float phi=pfid.phi_pim+180; //as per fid. cuts implemented by Isupov
				h3->Fill(pfid.theta_pim,phi);
			}
		}
		if (iSector==pfid.sector_p && pid.h10IdxP>0){
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			float phi=pfid.phi_p+180; //as per fid. cuts implemented by Isupov
			h1->Fill(pfid.theta_p,phi);
		}
		if (iSector==pfid.sector_pip && pid.h10IdxPip>0){
			TH2* h2 = (TH2*)hists[iSector]->At(1);
			float phi=pfid.phi_pip+180; //as per fid. cuts implemented by Isupov
			h2->Fill(pfid.theta_pip,phi);
		}
		if (iSector==pfid.sector_pim && pid.h10IdxPim>0 ){
			TH2* h3 = (TH2*)hists[iSector]->At(2);
			float phi=pfid.phi_pim+180; //as per fid. cuts implemented by Isupov
			h3->Fill(pfid.theta_pim,phi);
		}
	}
}

void DataAna::fillHistsPFidElast(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if (iSector==0) {
			if (pid_elast.h10IdxP>0){
				TH2* h1 = (TH2*)hists[iSector]->At(0);
				float phi=pfid_elast.phi_p+180; //as per fid. cuts implemented by Isupov
				h1->Fill(pfid_elast.theta_p,phi);
			}
		}
		if (iSector==pfid_elast.sector_p && pid_elast.h10IdxP>0){
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			float phi=pfid_elast.phi_p+180; //as per fid. cuts implemented by Isupov
			h1->Fill(pfid_elast.theta_p,phi);
		}
	}
}

void DataAna::fillHistsPFidElastNew(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if (iSector==0) {
			if (pid_elastnew.h10IdxP>0){
				TH2* h1 = (TH2*)hists[iSector]->At(0);
				float phi=pfid_elast.phi_p+180; //as per fid. cuts implemented by Isupov
				h1->Fill(pfid_elast.theta_p,phi);
			}
		}
		if (iSector==pfid_elast.sector_p && pid_elastnew.h10IdxP>0){
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			float phi=pfid_elast.phi_p+180; //as per fid. cuts implemented by Isupov
			h1->Fill(pfid_elast.theta_p,phi);
		}
	}
}

void DataAna::fillHistsPFidNew(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if (iSector==0) {
			if (pidnew.h10IdxP>0){
				TH2* h1 = (TH2*)hists[iSector]->At(0);
				float phi=pfid.phi_p+180; //as per fid. cuts implemented by Isupov
				h1->Fill(pfid.theta_p,phi);
			}
			if (pidnew.h10IdxPip>0){
				TH2* h2 = (TH2*)hists[iSector]->At(1);
				float phi=pfid.phi_pip+180; //as per fid. cuts implemented by Isupov
				h2->Fill(pfid.theta_pip,phi);
				
			}
			if (pidnew.h10IdxPim>0 ){
				TH2* h3 = (TH2*)hists[iSector]->At(2);
				float phi=pfid.phi_pim+180; //as per fid. cuts implemented by Isupov
				h3->Fill(pfid.theta_pim,phi);
			}
		}
		if (iSector==pfid.sector_p && pidnew.h10IdxP>0){
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			float phi=pfid.phi_p+180; //as per fid. cuts implemented by Isupov
			h1->Fill(pfid.theta_p,phi);
		}
		if (iSector==pfid.sector_pip && pidnew.h10IdxPip>0){
			TH2* h2 = (TH2*)hists[iSector]->At(1);
			float phi=pfid.phi_pip+180; //as per fid. cuts implemented by Isupov
			h2->Fill(pfid.theta_pip,phi);
		}
		if (iSector==pfid.sector_pim && pidnew.h10IdxPim>0 ){
			TH2* h3 = (TH2*)hists[iSector]->At(2);
			float phi=pfid.phi_pim+180; //as per fid. cuts implemented by Isupov
			h3->Fill(pfid.theta_pim,phi);
		}
	}
}

void DataAna::fillHistsEEff(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		/*if ( (iSector == 0) || (iSector == eeff.sector) ) {
		
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			h1->Fill(eeff.p, eeff.theta);
			
		}*/
		if (iSector==0) {
		
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			h1->Fill(eeff.p, eeff.theta);
			
		}
		if (iSector==eeff.sector) {
		
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			h1->Fill(eeff.p, eeff.theta);
			
		}
	}
}

void DataAna::fillHistsPEff(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if (iSector==0) {
			if (pid.h10IdxP>0){
				TH2* h1 = (TH2*)hists[iSector]->At(0);
				h1->Fill(peff.p_p, peff.theta_p);
			}
			if (pid.h10IdxPip>0){
				TH2* h2 = (TH2*)hists[iSector]->At(1);
				h2->Fill(peff.p_pip, peff.theta_pip);
				
			}
			if (pid.h10IdxPim>0 ){
				TH2* h3 = (TH2*)hists[iSector]->At(2);
				h3->Fill(peff.p_pim, peff.theta_pim);
			}
		}
		if (iSector==peff.sector_p && pid.h10IdxP>0){
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			h1->Fill(peff.p_p, peff.theta_p);
		}
		if (iSector==peff.sector_pip && pid.h10IdxPip>0){
			TH2* h2 = (TH2*)hists[iSector]->At(1);
			h2->Fill(peff.p_pip, peff.theta_pip);
		}
		if (iSector==peff.sector_pim && pid.h10IdxPim>0 ){
			TH2* h3 = (TH2*)hists[iSector]->At(2);
			h3->Fill(peff.p_pim, peff.theta_pim);
		}
	}
}

void DataAna::fillHistsPEffNew(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for(Int_t iSector=0;iSector<NSECTORS;iSector++){
		if (iSector==0) {
			if (pidnew.h10IdxP>0){
				TH2* h1 = (TH2*)hists[iSector]->At(0);
				h1->Fill(peff.p_p, peff.theta_p);
			}
			if (pidnew.h10IdxPip>0){
				TH2* h2 = (TH2*)hists[iSector]->At(1);
				h2->Fill(peff.p_pip, peff.theta_pip);
				
			}
			if (pidnew.h10IdxPim>0 ){
				TH2* h3 = (TH2*)hists[iSector]->At(2);
				h3->Fill(peff.p_pim, peff.theta_pim);
			}
		}
		if (iSector==peff.sector_p && pidnew.h10IdxP>0){
			TH2* h1 = (TH2*)hists[iSector]->At(0);
			h1->Fill(peff.p_p, peff.theta_p);
		}
		if (iSector==peff.sector_pip && pidnew.h10IdxPip>0){
			TH2* h2 = (TH2*)hists[iSector]->At(1);
			h2->Fill(peff.p_pip, peff.theta_pip);
		}
		if (iSector==peff.sector_pim && pidnew.h10IdxPim>0 ){
			TH2* h3 = (TH2*)hists[iSector]->At(2);
			h3->Fill(peff.p_pim, peff.theta_pim);
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

void DataAna::fillHistsPidMon(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for (int itrk=0;itrk<pidnew.ntrk;itrk++){//! loop over ntrk in the event
		int hIdx=-1;
		//! First fill PID hists for +ve trks
		if (pidnew.q[itrk]==1){//! Select +ve trks
			hIdx=0;
			for (int isctr=0;isctr<NSECTORS;isctr++){//! Loop over all sectors
				if ( (isctr == 0) || (isctr == pidnew.sector[itrk]-1) ) {//! Select appropriate sectors
					// Fill Beta v. p for all +ve trks
					TH2* h1 = (TH2*)hists[isctr]->At(hIdx);
					h1->Fill(pidnew.p[itrk], pidnew.b[itrk]);
					//! Fill dt v.p historam under proton hypothesis
					TH2* h2 = (TH2*)hists[isctr]->At(hIdx+1);
					h2->Fill(pidnew.p[itrk], pidnew.dt_p[itrk]);
					//! Now fill dt v.p historam under proton hypothesis
					TH2* h3 = (TH2*)hists[isctr]->At(hIdx+2);
					h3->Fill(pidnew.p[itrk], pidnew.dt_pip[itrk]);
				}
			}

		}
		//! NOw fill PID hists for -ve trks
		if (pidnew.q[itrk]==-1){//! Select -ve trks
			hIdx=3;
			for (int isctr=0;isctr<NSECTORS;isctr++){//! Loop over all sectors
				if ( (isctr == 0) || (isctr == pidnew.sector[itrk]-1) ) {//! Select appropriate sectors
					// Fill Beta v. p for all +ve trks
					TH2* h1 = (TH2*)hists[isctr]->At(hIdx);
					h1->Fill(pidnew.p[itrk], pidnew.b[itrk]);
					//! Fill dt v.p historam under pim hypothesis
					TH2* h2 = (TH2*)hists[isctr]->At(hIdx+1);
					h2->Fill(pidnew.p[itrk], pidnew.dt_pim[itrk]);
				}
			}

		}
	}
}

void DataAna::fillHistsPidCut(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for (int itrk=0;itrk<pidnew.ntrk;itrk++){//! loop over ntrk in the event
		int hIdx=-1;
		//! First fill PID hists for +ve trks
		if (pidnew.q[itrk]==1){//! Select +ve trks
			for (int isctr=0;isctr<NSECTORS;isctr++){//! Loop over all sectors
				if ( (isctr == 0) || (isctr == pidnew.sector[itrk]-1) ) {//! Select appropriate sectors
					// Fill PID hists for idtfd p
					if (pidnew.id[itrk]==PROTON && pidnew.h10IdxP>0){
						hIdx=0;
						TH2* h1 = (TH2*)hists[isctr]->At(hIdx);
						h1->Fill(pidnew.p[itrk], pidnew.b[itrk]);
						//! Fill dt v.p historam under proton hypothesis
						TH2* h2 = (TH2*)hists[isctr]->At(hIdx+1);
						h2->Fill(pidnew.p[itrk], pidnew.dt_p[itrk]);
					}
					// Fill PID hists for idtfd pip
					if (pidnew.id[itrk]==PIP && pidnew.h10IdxPip>0){
						hIdx=2;
						TH2* h1 = (TH2*)hists[isctr]->At(hIdx);
						h1->Fill(pidnew.p[itrk], pidnew.b[itrk]);
						//! Fill dt v.p historam under proton hypothesis
						TH2* h2 = (TH2*)hists[isctr]->At(hIdx+1);
						h2->Fill(pidnew.p[itrk], pidnew.dt_pip[itrk]);
					}
				}
			}
		}
		//! NOw fill PID hists for -ve trks
		if (pidnew.q[itrk]==-1){//! Select -ve trks
			hIdx=3;
			for (int isctr=0;isctr<NSECTORS;isctr++){//! Loop over all sectors
				if ( (isctr == 0) || (isctr == pidnew.sector[itrk]-1) ) {//! Select appropriate sectors
					// Fill PID hists for idtfd pim
					if (pidnew.id[itrk]==PIM && pidnew.h10IdxPim>0){
						hIdx=4;
						TH2* h1 = (TH2*)hists[isctr]->At(hIdx);
						h1->Fill(pidnew.p[itrk], pidnew.b[itrk]);
						//! Fill dt v.p historam under proton hypothesis
						TH2* h2 = (TH2*)hists[isctr]->At(hIdx+1);
						h2->Fill(pidnew.p[itrk], pidnew.dt_pim[itrk]);
					}
				}
			}
		}
	}
}

void DataAna::fillHistsPidElast(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	Int_t hIdx = 0;
	if(pid_elast.h10IdxP > 0){
		hIdx = 0;
		for(Int_t iSector=0;iSector<NSECTORS;iSector++){
				if ( (iSector == 0) || (iSector == pid_elast.sectorP) ) {
			
				TH2* h1 = (TH2*)hists[iSector]->At(hIdx);
				h1->Fill(pid_elast.pP, pid_elast.betaP);
				TH2* h2 = (TH2*)hists[iSector]->At(hIdx+1);
				h2->Fill(pid_elast.pP, pid_elast.dtP);
		
				TH2* h3 = (TH2*)hists[iSector]->At(hIdx+2);
				h3->Fill(pid_elast.P_ec_ei, pid_elast.P_ec_eo);
				TH2* h4 = (TH2*)hists[iSector]->At(hIdx+3);
				h4->Fill(pid_elast.P_ec_ei/pid_elast.pP, pid_elast.P_ec_eo/pid_elast.pP);
				TH2* h5 = (TH2*)hists[iSector]->At(hIdx+4);
				h5->Fill(pid_elast.pP, pid_elast.P_etot/pid_elast.pP);
				TH1* h6 = (TH1*)hists[iSector]->At(hIdx+5);
				h6->Fill(pid_elast.P_nphe);
				TH2* h7 = (TH2*)hists[iSector]->At(hIdx+6);
				h7->Fill(pid_elast.P_cc_segm, pid_elast.P_cc_theta);
			}
		}
	}
}

void DataAna::fillHistsPidElastMon(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for (int itrk=0;itrk<pid_elastnew.ntrk;itrk++){//! loop over ntrk in the event
		int hIdx=-1;
		//! Fill PID hists for +ve trks (there should only be +ve tracks!)
		if (pid_elastnew.q[itrk]==1){//! Select +ve trks (there should only be +ve tracks!)
			hIdx=0;
			for (int isctr=0;isctr<NSECTORS;isctr++){//! Loop over all sectors
				if ( (isctr == 0) || (isctr == pid_elastnew.sector[itrk]-1) ) {//! Select appropriate sectors
					// Fill Beta v. p for all +ve trks
					TH2* h1 = (TH2*)hists[isctr]->At(hIdx);
					h1->Fill(pid_elastnew.p[itrk], pid_elastnew.b[itrk]);
					//! Fill dt v.p historam under proton hypothesis
					TH2* h2 = (TH2*)hists[isctr]->At(hIdx+1);
					h2->Fill(pid_elastnew.p[itrk], pid_elastnew.dt_p[itrk]);
					//! Now fill dt v.p historam under proton hypothesis
					TH2* h3 = (TH2*)hists[isctr]->At(hIdx+2);
					h3->Fill(pid_elastnew.p[itrk], pid_elastnew.dt_pip[itrk]);
				}
			}

		}
	}
}

void DataAna::fillHistsPidElastCut(TObjArray** hists, Bool_t useMc /* = kFALSE */)
{
	for (int itrk=0;itrk<pid_elastnew.ntrk;itrk++){//! loop over ntrk in the event
		int hIdx=-1;
		//! First fill PID hists for +ve trks
		if (pid_elastnew.q[itrk]==1){//! Select +ve trks
			for (int isctr=0;isctr<NSECTORS;isctr++){//! Loop over all sectors
				if ( (isctr == 0) || (isctr == pid_elastnew.sector[itrk]-1) ) {//! Select appropriate sectors
					// Fill PID hists for idtfd p
					if (pid_elastnew.id[itrk]==PROTON && pid_elastnew.h10IdxP>0){
						hIdx=0;
						TH2* h1 = (TH2*)hists[isctr]->At(hIdx);
						h1->Fill(pid_elastnew.p[itrk], pid_elastnew.b[itrk]);
						//! Fill dt v.p historam under proton hypothesis
						TH2* h2 = (TH2*)hists[isctr]->At(hIdx+1);
						h2->Fill(pid_elastnew.p[itrk], pid_elastnew.dt_p[itrk]);
					}
				}
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
		hIdx=0;
		TH1* h1 = (TH1*)hists->At(hIdx);
		h1->Fill(tp->mm2ppippim);
		TH1* h2 = (TH1*)hists->At(hIdx+1);
		h2->Fill(tp->mmppippim);

		hIdx=8;
		TH2* h3 = (TH2*)hists->At(hIdx);
		h3->Fill(tp->W,tp->mm2ppippim);
		TH2* h4 = (TH2*)hists->At(hIdx+1);
		h4->Fill(tp->W,tp->mmppippim);
	}
	
    if ((useMc || (h10idxP>0 && h10idxPip>0 && h10idxPim==-1)) ) {
		hIdx=2;
		TH1* h1 = (TH1*)hists->At(hIdx);
		h1->Fill(tp->mm2ppip);
		TH1* h2 = (TH1*)hists->At(hIdx+1);
		h2->Fill(tp->mmppip);

		hIdx=10;
		TH2* h3 = (TH2*)hists->At(hIdx);
		h3->Fill(tp->W,tp->mm2ppip);
		TH2* h4 = (TH2*)hists->At(hIdx+1);
		h4->Fill(tp->W,tp->mmppip);
	}
	
	if ((useMc || (h10idxP>0 && h10idxPim>0 && h10idxPip==-1)) ) {
		hIdx=4;
		TH1* h1 = (TH1*)hists->At(hIdx);
		h1->Fill(tp->mm2ppim);
		TH1* h2 = (TH1*)hists->At(hIdx+1);
		h2->Fill(tp->mmppim);

		hIdx=12;
		TH2* h3 = (TH2*)hists->At(hIdx);
		h3->Fill(tp->W,tp->mm2ppim);
		TH2* h4 = (TH2*)hists->At(hIdx+1);
		h4->Fill(tp->W,tp->mmppim);
	}
	
	if ((useMc || (h10idxPip>0 && h10idxPim>0 && h10idxP==-1)) ) {
		hIdx = 6;
		TH1* h1 = (TH1*)hists->At(hIdx);
		h1->Fill(tp->mm2pippim);
		TH1* h2 = (TH1*)hists->At(hIdx+1);
		h2->Fill(tp->mmpippim);

		hIdx=14;
		TH2* h3 = (TH2*)hists->At(hIdx);
		h3->Fill(tp->W,tp->mm2pippim);
		TH2* h4 = (TH2*)hists->At(hIdx+1);
		h4->Fill(tp->W,tp->mmpippim);
	}
	
}

void DataAna::fillHistsMMElastic(TObjArray *hists, Bool_t useMc /* = kFALSE */)
{
	DataElastic *de = &dElast;
	if (useMc) de = &dElast_ST;
		
	TH2* h1 = (TH2*)hists->At(0);
	h1->Fill(de->W,de->MMep);
	TH1* h2 = (TH1*)hists->At(1);
	h2->Fill(de->W);
}

void DataAna::fillYields(TObjArray **hists, Float_t w, Bool_t useMc /* = kFALSE */)
{
	//! Get index of Crw-W bin
	/*int crswidx=GetCrsWBinIdx(2.950);
	Info("DataAna::fillYields()","crswidx for %f = %d",2.950,crswidx);*/
	int iw=GetCrsWBinIdx(w);
	/*Info("DataAna::fillYields()","crswidx for %f = %d",w,iw);*/
	if (iw==9999) return;

	//! Get h8[w]
	TObjArray* h8=hists[iw];

	Data2pi *tp = &d2pi;
	if (useMc) tp = &d2pi_mc;

	//cout<<"alphas2="<<tp->alpha_1<<":"<<tp->alpha_2<<":"<<tp->alpha_3<<endl;
	THnSparse* h8_1 = (THnSparse*)h8->At(0);
	//Double_t coord1[] = { tp->h, tp->Q2, tp->W, tp->varset1.M1, tp->varset1.M2, tp->varset1.theta, tp->varset1.phi, tp->varset1.alpha };
	Double_t coord1[] = { tp->h,tp->Q2,tp->W,tp->M_ppip,tp->M_pippim,tp->theta_cms_pim,tp->phi_cms_pim,tp->alpha_1 };
	h8_1->Fill(coord1);
	
	THnSparse* h8_2 = (THnSparse*)h8->At(1);
	//Double_t coord2[] = { tp->h, tp->Q2, tp->W, tp->varset2.M1, tp->varset2.M2, tp->varset2.theta, tp->varset2.phi, tp->varset2.alpha  };
	Double_t coord2[] = { tp->h,tp->Q2,tp->W,tp->M_ppip,tp->M_pippim,tp->theta_cms_p,tp->phi_cms_p,tp->alpha_2};
	h8_2->Fill(coord2);
	
	THnSparse* h8_3 = (THnSparse*)h8->At(2);
	//Double_t coord3[] = { tp->h, tp->Q2, tp->W, tp->varset3.M1, tp->varset3.M2, tp->varset3.theta, tp->varset3.phi, tp->varset3.alpha  };
	Double_t coord3[] = { tp->h,tp->Q2,tp->W,tp->M_ppip,tp->M_ppim,tp->theta_cms_pip,tp->phi_cms_pip,tp->alpha_3};
	h8_3->Fill(coord3);

	return;
}

void DataAna::fillYieldsElastic(TObjArray *hists, Bool_t useMc /* = kFALSE */)
{
	DataElastic *de = &dElast;
	if (useMc) de = &dElast_ST;

	THnSparse* h2 = (THnSparse*)hists->At(0);
	Double_t coord[] = {de->theta_e,de->phi_e};
	h2->Fill(coord);
}

void DataAna::addBranches_Data2pi(TTree* t, Bool_t useMc/*=kFALSE*/){
	Data2pi *tp = &d2pi;
	if (useMc) tp = &d2pi_mc;

	//! Initial Beam Energy
	t->Branch("p_e0",&tp->p_e0);
	//! Reconstructed Kinematics 
	//! for e',p',p,pip,pim at e' vertex
	t->Branch("p_e",&tp->p_e);
	t->Branch("p_p",&tp->p_p);
	t->Branch("p_pip",&tp->p_pip);
	t->Branch("p_pim",&tp->p_pim);
	t->Branch("theta_e",&tp->theta_e);
	t->Branch("theta_p",&tp->theta_p);
	t->Branch("theta_pip",&tp->theta_pip);
	t->Branch("theta_pim",&tp->theta_pim);
	t->Branch("phi_e",&tp->phi_e);
	t->Branch("phi_p",&tp->phi_p);
	t->Branch("phi_pip",&tp->phi_pip);
	t->Branch("phi_pim",&tp->phi_pim);
	//! Reconstructed e' Vertex
	t->Branch("vx_e",&tp->vx_e);
	t->Branch("vx_p",&tp->vx_p);
	t->Branch("vx_pip",&tp->vx_pip);
	t->Branch("vx_pim",&tp->vx_pim);
	t->Branch("vy_e",&tp->vy_e);
	t->Branch("vy_p",&tp->vy_p);
	t->Branch("vy_pip",&tp->vy_pip);
	t->Branch("vy_pim",&tp->vy_pim);
	t->Branch("vz_e",&tp->vz_e);
	t->Branch("vz_p",&tp->vz_p);
	t->Branch("vz_pip",&tp->vz_pip);
	t->Branch("vz_pim",&tp->vz_pim);
	//! Q2, W
	t->Branch("Q2",&tp->Q2);
	t->Branch("W",&tp->W);
	//! Helicity
	t->Branch("h",&tp->h);
	//! {MMs}
	t->Branch("mm2ppippim",&tp->mm2ppippim);
	t->Branch("mmppippim",&tp->mmppippim);
	t->Branch("mm2ppip",&tp->mm2ppip);
	t->Branch("mmppip",&tp->mmppip);
	t->Branch("mm2ppim",&tp->mm2ppim);
	t->Branch("mmppim",&tp->mmppim);
	t->Branch("mm2pippim",&tp->mm2pippim);
	t->Branch("mmpippim",&tp->mmpippim);
	//! Topology
	t->Branch("top",&tp->top);
	//! Varsets
	t->Branch("M_ppip",&tp->M_ppip);
	t->Branch("M_ppim",&tp->M_ppim);
	t->Branch("M_pippim",&tp->M_pippim);
	t->Branch("theta_cms_p",&tp->theta_cms_p);
	t->Branch("theta_cms_pip",&tp->theta_cms_pip);
	t->Branch("theta_cms_pim",&tp->theta_cms_pim);
	t->Branch("phi_cms_p",&tp->phi_cms_p);
	t->Branch("phi_cms_pip",&tp->phi_cms_pip);
	t->Branch("phi_cms_pim",&tp->phi_cms_pim);
	t->Branch("alpha_1",&tp->alpha_1);
	t->Branch("alpha_2",&tp->alpha_2);
	t->Branch("alpha_3",&tp->alpha_3);
}

void DataAna::addBranches_DataElastic(TTree* t, Bool_t useMc/*=kFALSE*/){
	DataElastic* de=&dElast;
	if (useMc) de=&dElast_ST;
	//! gpart
	t->Branch("gpart",&de->gpart);
	//! Q2, W
	t->Branch("Q2",&de->Q2);
	t->Branch("W", &de->W);
	//! MMs for Event Selection
	t->Branch("MMep",  &de->MMep);
	t->Branch("MM2ep", &de->MM2ep);
	//! Reconstructed Kinematics 
	//! for e',p',p,pip,pim at e' vertex
	t->Branch("p_e",&de->p_e);
	t->Branch("p_p",&de->p_p);
	t->Branch("theta_e",&de->theta_e);
	t->Branch("theta_p",&de->theta_p);
	t->Branch("phi_e",&de->phi_e);
	t->Branch("phi_p",&de->phi_p);
	//! Reconstructed e' Vertex
	t->Branch("vx_e",&de->vx_e);
	t->Branch("vx_p",&de->vx_p);
	t->Branch("vy_e",&de->vy_e);
	t->Branch("vy_p",&de->vy_p);
	t->Branch("vz_e",&de->vz_e);
	t->Branch("vz_p",&de->vz_p);
	//! EID
	//t->Branch("nphe",&de->nphe);
}

void DataAna::addBranches_DataEid(TTree* t){
	t->Branch("sector_e",&eid.sector);
    t->Branch("beta_e",&eid.beta);
    t->Branch("betaStr_e",&eid.betaStrE);
	t->Branch("dt_e",&eid.dtE);
	t->Branch("p_e_eid",&eid.p); //so that it does not coflict with p_e from DataElast 
	//from EC
	t->Branch("ec_ei",&eid.ec_ei);
	t->Branch("ec_eo",&eid.ec_eo);
	t->Branch("etot",&eid.etot);
	//from CC
	t->Branch("nphe",&eid.nphe);
	t->Branch("cc_segm",&eid.cc_segm);
	t->Branch("cc_theta",&eid.cc_theta);
}

void DataAna::addBranches_DataPid(TTree* t){
	t->Branch("sector_p",&pid.sectorP);
    t->Branch("sector_pip",&pid.sectorPip);
    t->Branch("sector_pim",&pid.sectorPim);
    t->Branch("beta_p",&pid.betaP);
    t->Branch("beta_pip",&pid.betaPip);
    t->Branch("beta_pim",&pid.betaPim);
    t->Branch("betaStr_p",&pid.betaStrP);
    t->Branch("betaStr_pip",&pid.betaStrPip);
    t->Branch("betaStr_pim",&pid.betaStrPim);
    t->Branch("dt_p",&pid.dtP);
    t->Branch("dt_pip",&pid.dtPip);
    t->Branch("dt_pim",&pid.dtPim);
    t->Branch("p_p_pid",&pid.pP);
    t->Branch("p_pip_pid",&pid.pPip);
    t->Branch("p_pim_pid",&pid.pPim);
    /*CC, EC signature*/
    //from EC
    t->Branch("ec_ei_p",&pid.P_ec_ei);
    t->Branch("ec_ei_pip",&pid.Pip_ec_ei);
    t->Branch("ec_ei_pim",&pid.Pim_ec_ei);
    t->Branch("ec_eo_p",&pid.P_ec_eo);
    t->Branch("ec_eo_pip",&pid.Pip_ec_eo);
    t->Branch("ec_eo_pim",&pid.Pim_ec_eo);
    t->Branch("etot_p",&pid.P_etot);
    t->Branch("etot_pip",&pid.Pip_etot);
    t->Branch("etot_pim",&pid.Pim_etot);
    //from CC
    t->Branch("nphe_p",&pid.P_nphe);
    t->Branch("nphe_pip",&pid.Pip_nphe);
    t->Branch("nphe_pim",&pid.Pim_nphe);
    t->Branch("cc_segm_p",&pid.P_cc_segm);
    t->Branch("cc_segm_pip",&pid.Pip_cc_segm);
    t->Branch("cc_segm_pim",&pid.Pim_cc_segm);
    t->Branch("cc_theta_p",&pid.P_cc_theta);
    t->Branch("cc_theta_pip",&pid.Pip_cc_theta);
    t->Branch("cc_theta_pim",&pid.Pim_cc_theta);
}

void DataAna::addBranches_DataPidNew(TTree* t){
	t->Branch("l_e",&pidnew.l_e,"l_e/F");
	t->Branch("t_e",&pidnew.t_e,"t_e/F");
	t->Branch("t_off",&pidnew.t_off,"t_off/F");
	t->Branch("ntrk",&pidnew.ntrk,"ntrk/I");
	t->Branch("q",pidnew.q,"q[ntrk]/I");
	t->Branch("dc",pidnew.dc,"dc[ntrk]/I");
	t->Branch("sc",pidnew.sc,"sc[ntrk]/I");
	t->Branch("q",pidnew.q,"q[ntrk]/I");
	t->Branch("p",pidnew.p,"p[ntrk]/F");
	t->Branch("l",pidnew.l,"l[ntrk]/F");
	t->Branch("t",pidnew.t,"t[ntrk]/F");
	t->Branch("b",pidnew.b,"b[ntrk]/F");
	t->Branch("id",pidnew.id,"id[ntrk]/I");
	t->Branch("sector",pidnew.sector,"sector[ntrk]/I");
	t->Branch("h10_idx",pidnew.h10_idx,"h10_idx[ntrk]/I");
	t->Branch("b_p",pidnew.b_p,"b_p[ntrk]/F");
	t->Branch("b_pip",pidnew.b_pip,"b_pip[ntrk]/F");
	t->Branch("b_pim",pidnew.b_pim,"b_pim[ntrk]/F");
	t->Branch("dt_p",pidnew.dt_p,"dt_p[ntrk]/F");
	t->Branch("dt_pip",pidnew.dt_pip,"dt_pip[ntrk]/F");
	t->Branch("dt_pim",pidnew.dt_pim,"dt_pim[ntrk]/F");
}

void DataAna::addBranches_DataPidElast(TTree* t){
	t->Branch("sector_p",&pid_elast.sectorP);
    t->Branch("beta_p",&pid_elast.betaP);
    t->Branch("betaStr_p",&pid_elast.betaStrP);
    t->Branch("dt_p",&pid_elast.dtP);
    t->Branch("p_p_pid",&pid_elast.pP);
    /*CC, EC signature*/
    //from EC
    t->Branch("ec_ei_p",&pid_elast.P_ec_ei);
    t->Branch("ec_eo_p",&pid_elast.P_ec_eo);
    t->Branch("etot_p",&pid_elast.P_etot);
    //from CC
    t->Branch("nphe_p",&pid_elast.P_nphe);
    t->Branch("cc_segm_p",&pid_elast.P_cc_segm);
    t->Branch("cc_theta_p",&pid_elast.P_cc_theta);
}