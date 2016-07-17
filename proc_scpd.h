#ifndef PROCSCPD_H
#define PROCSCPD_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "wrpr_cut_sc_pd_e16.h" //! [12-11-15] EI's sc_pd cuts
#include <TMath.h>

using namespace TMath;

/*
 [07-13-16]
 + Added to emulate E16's lite:scpd cuts and therefore lite code may be reproduced here in varying degrees.
	 
 + Processor must follow ProcEid:ProcPidNew:ProcMom

 + Since 'lite:eff' is implemented only for t2'(e,p,pip,pim-ignore), the fully exclusive selection of the hvy machinery
   (t1,t2,t3 and t4) is broken down into the following 2 types of selection, whose criterion if satisfied, the events are passed.
   NOTE that only in the 1st selection, i.e. t2', are the 'lite:eff' cuts valid and therefore applied:
    1. t2'(~t1+t2) criterion: 
       If 
       + dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0
        OR
       + dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0
        AND
       + each of e,p& pip are in efficient region of detector
 	
	2. t3+t4: If in these topologies, pass events without requiring 'lite:eff' cuts since these do not exist for pim:
	  If 
	  + dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPim>0 
	    OR
	  + dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0 

*/
class ProcScpd : public EpProcessor {

public:
	ProcScpd(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		     Bool_t monitor=kFALSE,Bool_t monitorOnly=kFALSE, Bool_t use_scpd_at_mod=kFALSE);
	~ProcScpd();
	
	void handle();
	//void write();
	
protected:	
	
	//! for hevtsum
	static const int NUM_EVTCUTS=5;
	enum {EVT_NULL, EVT_TOT, EVT_E_PASS, EVT_P_PASS, EVT_PIP_PASS, EVT_E_AND_P_AND_PIP_PASS};
	
	//! _hscpd[3][2] for e,p,pip:prec,pstc
	TString* _scpd_prtcl_names;
	TH1D*** _hpdl;
	//! _hthetavp[3][6][2] for e,p,pip:sector:prec,pstc
	TH2D**** _hthetavp2;

	bool _use_scpd_at_mod;

	void updateScpd();
	bool pass_scpd();
	bool is_scpd_bad(TString prtcl_name);
	int get_sector(float phi);
};

ProcScpd::ProcScpd(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
                   Bool_t monitor/* = kFALSE*/,Bool_t monitorOnly /*= kFALSE*/,
				   Bool_t use_scpd_at_mod/* = kFALSE*/)
                   :EpProcessor(td, dataH10, dataAna, monitor, monitorOnly)
{
	_use_scpd_at_mod=use_scpd_at_mod;
	if (_use_scpd_at_mod){
		Info("ProcScpd::ProcScpd()", "Will use_scpd_at_mod");
	}

	td->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT_TOT,"total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_E_PASS,"e^{-} inscpd");
	hevtsum->GetXaxis()->SetBinLabel(EVT_P_PASS,"p inscpd");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIP_PASS,"#pi^{+} inscpd");
	hevtsum->GetXaxis()->SetBinLabel(EVT_E_AND_P_AND_PIP_PASS,"e+p+#pi^{+} inscpd");
	hevtsum->SetMinimum(0.);

	//! cut monitor histograms
	//! _hpdl[3][2] for e,p,pip:prec,pstc
	_scpd_prtcl_names=new TString[3];
	_scpd_prtcl_names[0]="e";
	_scpd_prtcl_names[1]="p";
	_scpd_prtcl_names[2]="pip";
	_hpdl=new TH1D**[3]; //!e,p,pip
	for (int i=0;i<3;i++){
		_hpdl[i]=new TH1D*[2]; //!prec,pstc
		for (int j=0;j<2;j++){
			TString name_cut_lvl;
			if      (j==0) name_cut_lvl="prec";
			else if (j==1) name_cut_lvl="pstc";

			TString name=TString::Format("pdl_%s_%s",
				                         _scpd_prtcl_names[i].Data(),name_cut_lvl.Data());
			_hpdl[i][j]=new TH1D(name,name,700,0,700);
		}
	}
	//! _hthetavp2[3][6][2] for e,p,pip:sector:prec,pstc
	_hthetavp2=new TH2D***[3]; //!e,p,pip
	for (int i=0;i<3;i++){
		_hthetavp2[i]=new TH2D**[6]; //!sector
		for (int j=0;j<6;j++){
			_hthetavp2[i][j]=new TH2D*[2];//prec,pstc
			for (int k=0;k<2;k++){
				TString name_cut_lvl;
				if      (k==0) name_cut_lvl="prec";
				else if (k==1) name_cut_lvl="pstc";

				TString name=TString::Format("thetavp2_%s_s%d_%s",
					                         _scpd_prtcl_names[i].Data(),j+1,name_cut_lvl.Data());
				if      (_scpd_prtcl_names[i]=="e")   _hthetavp2[i][j][k]=new TH2D(name,name,100,0,5,100,0,60);
				else if (_scpd_prtcl_names[i]=="p")   _hthetavp2[i][j][k]=new TH2D(name,name,100,0,4,100,0,60);
				else if (_scpd_prtcl_names[i]=="pip") _hthetavp2[i][j][k]=new TH2D(name,name,100,0,3,100,0,120);	
			}
		}
	}
}

ProcScpd::~ProcScpd()
{
	delete[] _scpd_prtcl_names;
	delete[] _hpdl;
	delete[] _hthetavp2;
}

void ProcScpd::handle()
{
	//Info("ProcScpd::handle()", "");
	pass = kFALSE;
	hevtsum->Fill(EVT_TOT);

	updateScpd();

	bool in_scpd=kFALSE;
	if (dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0) { //! t1 and t2 pim does not matter; dAna->pidnew.h10IdxPim>0) {
		in_scpd=pass_scpd();
		if (in_scpd) {
			hevtsum->Fill(EVT_E_AND_P_AND_PIP_PASS);
		}
	}else if (dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip==0 && dAna->pidnew.h10IdxPim>0){//! t3
		in_scpd=kTRUE;
	}else if (dAna->pidnew.h10IdxP==0 && dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0){//! t4
		in_scpd=kTRUE;
	}
		
	if (in_scpd)
	{
		pass = kTRUE;
		EpProcessor::handle();
	}
}

/* Currently implemented for top2'*/
bool ProcScpd::pass_scpd(){
	//Info("h10looper_2pi::pass_scpd()","");
	bool ret=kFALSE;
	bool e_inscpd=kFALSE;
	bool p_inscpd=kFALSE;
	bool pip_inscpd=kFALSE;

	DataScpd* dscpd = &dAna->scpd;

	int sector_e =dscpd->sector_e;
	float p_e    =dscpd->p_e;
	float theta_e=dscpd->theta_e;
	int sc_pd_e  =dscpd->sc_pd_e;

	int sector_p =dscpd->sector_p;
	float p_p    =dscpd->p_p;
	float theta_p=dscpd->theta_p;
	int sc_pd_p  =dscpd->sc_pd_p;

	int sector_pip =dscpd->sector_pip;
	float p_pip    =dscpd->p_pip;
	float theta_pip=dscpd->theta_pip;
	int sc_pd_pip  =dscpd->sc_pd_pip;

	//Info("ProcScpd::pass_scpd()","sectors:%d, %d, %d",sector_e,sector_p,sector_pip);

	//!prec
	_hpdl[0][0]->Fill(sc_pd_e);
	_hpdl[1][0]->Fill(sc_pd_p);
	_hpdl[2][0]->Fill(sc_pd_pip);
	_hthetavp2[0][sector_e-1][0]->Fill(p_e,theta_e);
	_hthetavp2[1][sector_p-1][0]->Fill(p_p,theta_p);
	_hthetavp2[2][sector_pip-1][0]->Fill(p_pip,theta_pip);

	//! check if e- inscpd
	//Info("ProcScpd::pass_scpd()","sc_pd_e=%d",sc_pd_e);
	if (!is_scpd_bad("e")){
		hevtsum->Fill(EVT_E_PASS);
		e_inscpd=kTRUE;
		_hpdl[0][1]->Fill(sc_pd_e);
		_hthetavp2[0][sector_e-1][1]->Fill(p_e,theta_e);
	}
	//! check if p inscpd
	//Info("ProcScpd::pass_scpd()","sc_pd_p=%d",sc_pd_p);
	if (!is_scpd_bad("p")){
		hevtsum->Fill(EVT_P_PASS);
		p_inscpd=kTRUE;
		_hpdl[1][1]->Fill(sc_pd_p);
		_hthetavp2[1][sector_p-1][1]->Fill(p_p,theta_p);
	}
	//! check if pip inscpd
	//Info("ProcScpd::pass_scpd()","sc_pd_pip=%d",sc_pd_pip);
	if (!is_scpd_bad("pip")){
		hevtsum->Fill(EVT_PIP_PASS);
		pip_inscpd=kTRUE;
		_hpdl[2][1]->Fill(sc_pd_pip);
		_hthetavp2[2][sector_pip-1][1]->Fill(p_pip,theta_pip);
	}

	ret=e_inscpd && p_inscpd && pip_inscpd;
	return ret;
}

void ProcScpd::updateScpd() {
	TLorentzVector lvE,lvP,lvPip,lvPim;

	//! e
	Float_t mom = dH10->p[0];
	Float_t px = mom*dH10->cx[0];
	Float_t py = mom*dH10->cy[0];
	Float_t pz = mom*dH10->cz[0];
	lvE.SetPxPyPzE(px,py,pz,Sqrt(mom*mom+MASS_E*MASS_E));
	
	//! *** sector ***
	//dAna->eff.sector_e=dH10->sc_sect[dH10->sc[0]-1];	
	float phi_tmp_e=lvE.Phi()*TMath::RadToDeg();// [-180,180]
	float phi_e=phi_tmp_e<-30?phi_tmp_e+360:phi_tmp_e; // [-30,330]
	dAna->scpd.sector_e=get_sector(phi_e);
	//! ******	
	dAna->scpd.theta_e=lvE.Theta()*RadToDeg();
	dAna->scpd.p_e=lvE.P();
	dAna->scpd.sc_pd_e=100*dH10->dc_sect[dH10->dc[0]-1]+dH10->sc_pd[dH10->sc[0]-1];

	//!p
	if (dAna->pidnew.h10IdxP>0) {
		Double_t mom = dH10->p[dAna->pidnew.h10IdxP];
		Double_t px = mom*dH10->cx[dAna->pidnew.h10IdxP];
		Double_t py = mom*dH10->cy[dAna->pidnew.h10IdxP];
		Double_t pz = mom*dH10->cz[dAna->pidnew.h10IdxP];
		Double_t energy = Sqrt(mom*mom+MASS_P*MASS_P);
		lvP.SetPxPyPzE(px,py,pz,energy);

		//! *** sector ***
		//dAna->eff.sector_p=dH10->sc_sect[dH10->sc[dAna->pidnew.h10IdxP]-1];
		float phi_tmp_p=lvP.Phi()*TMath::RadToDeg();// [-180,180]
		float phi_p=phi_tmp_p<-30?phi_tmp_p+360:phi_tmp_p; // [-30,330]
		dAna->scpd.sector_p=get_sector(phi_p);
		//! ******
		dAna->scpd.theta_p=lvP.Theta()*RadToDeg();
		dAna->scpd.p_p=lvP.P();
		dAna->scpd.sc_pd_p=100*dH10->dc_sect[dH10->dc[dAna->pidnew.h10IdxP]-1]+dH10->sc_pd[dH10->sc[dAna->pidnew.h10IdxP]-1];
	}
	//!pip
	if (dAna->pidnew.h10IdxPip>0) {
		Double_t mom = dH10->p[dAna->pidnew.h10IdxPip];
		Double_t px = mom*dH10->cx[dAna->pidnew.h10IdxPip];
		Double_t py = mom*dH10->cy[dAna->pidnew.h10IdxPip];
		Double_t pz = mom*dH10->cz[dAna->pidnew.h10IdxPip];
		Double_t energy = Sqrt(mom*mom+MASS_PIP*MASS_PIP);
		lvPip.SetPxPyPzE(px,py,pz,energy);

		//! *** sector ***
		//dAna->eff.sector_pip=dH10->sc_sect[dH10->sc[dAna->pidnew.h10IdxPip]-1];
		float phi_tmp_pip=lvPip.Phi()*TMath::RadToDeg();// [-180,180]
		float phi_pip=phi_tmp_pip<-30?phi_tmp_pip+360:phi_tmp_pip; // [-30,330]
		dAna->scpd.sector_pip=get_sector(phi_pip);
		//! ******
		dAna->scpd.theta_pip=lvPip.Theta()*RadToDeg();
		dAna->scpd.p_pip=lvPip.P();
		dAna->scpd.sc_pd_pip=100*dH10->dc_sect[dH10->dc[dAna->pidnew.h10IdxPip]-1]+dH10->sc_pd[dH10->sc[dAna->pidnew.h10IdxPip]-1];
	}
}

bool ProcScpd::is_scpd_bad(TString prtcl_name){
	//Info("h10looper_e1f::pass_scpd()","");
	bool ret=kFALSE;

	DataScpd* dscpd = &dAna->scpd;

	int sc_pd_e  =dscpd->sc_pd_e;

	int sc_pd_p  =dscpd->sc_pd_p;

	int sc_pd_pip  =dscpd->sc_pd_pip;
	
	if (dH10->expt=="e1f"){
		ret=kTRUE;
	}else if (dH10->expt=="e16"){
		int at_mod=0;
		if (_use_scpd_at_mod) at_mod=1;

		if (prtcl_name=="e") {
			ret=is_scpd_bad_e16(sc_pd_e,at_mod);
		}else if (prtcl_name=="p") {
			ret=is_scpd_bad_e16(sc_pd_p,at_mod);
		}else if (prtcl_name=="pip") {
			ret=is_scpd_bad_e16(sc_pd_pip,at_mod);
		}
	}
	return ret;
}

/*
phi has to be [-30,330] for the following function to work.
*/
int ProcScpd::get_sector(float phi){
	int sector=-1;
	if ( (phi>=PHI_FULL_SECTOR[0][0]  && phi<PHI_FULL_SECTOR[0][1]) || phi==PHI_FULL_SECTOR[5][1]){
		sector=1;
	}else if (phi>=PHI_FULL_SECTOR[1][0] && phi<PHI_FULL_SECTOR[1][1]){
		sector=2;
	}else if (phi>=PHI_FULL_SECTOR[2][0] && phi<PHI_FULL_SECTOR[2][1]){
		sector=3;
	}else if (phi>=PHI_FULL_SECTOR[3][0] && phi<PHI_FULL_SECTOR[3][1]){
		sector=4;
	}else if (phi>=PHI_FULL_SECTOR[4][0] && phi<PHI_FULL_SECTOR[4][1]){
		sector=5;
	}else if (phi>=PHI_FULL_SECTOR[5][0] && phi<PHI_FULL_SECTOR[5][1]){
		sector=6;
	}else{
		Info("ProcScpd::get_sector","sector not found for phi=%.2f. Returning sector=-1",phi);
		//Info("ProcEff::get_sector","file=%s\n",fChain->GetCurrentFile()->GetName());
	}
	return sector;
}

#endif // PROCSCPD_H
