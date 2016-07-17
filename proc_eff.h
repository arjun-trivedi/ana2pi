#ifndef PROCEFF_H
#define PROCEFF_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "wrpr_cut_theta_vs_p_e16.h" //! [12-10-15] EI's theta_vs_p cuts
#include <TMath.h>

using namespace TMath;

/*
 [07-13-16]
 + Added to emulate E16's lite:eff cuts and therefore lite code may be reproduced here in varying degrees.
	 
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
class ProcEff : public EpProcessor {

public:
	ProcEff(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		     Bool_t monitor=kFALSE,Bool_t monitorOnly=kFALSE, Bool_t use_eff_at_mod=kFALSE);
	~ProcEff();
	
	void handle();
	//void write();
	
protected:	
	
	//! for hevtsum
	static const int NUM_EVTCUTS=5;
	enum {EVT_NULL, EVT_TOT, EVT_E_PASS, EVT_P_PASS, EVT_PIP_PASS, EVT_E_AND_P_AND_PIP_PASS};
	
	//! _hthetavp[3][6][2] for e,p,pip:sector:prec,pstc
 	TString* _eff_prtcl_names;
	TH2D**** _hthetavp;

	bool _use_eff_at_mod;

	void updateEff();
	bool pass_eff();
	bool pass_theta_vs_p(TString prtcl_name);
	int get_sector(float phi);
};

ProcEff::ProcEff(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
                   Bool_t monitor/* = kFALSE*/,Bool_t monitorOnly /*= kFALSE*/,
				   Bool_t use_eff_at_mod/* = kFALSE*/)
                   :EpProcessor(td, dataH10, dataAna, monitor, monitorOnly)
{
	_use_eff_at_mod=use_eff_at_mod;
	if (_use_eff_at_mod){
		Info("ProcEff::ProcEff()", "Will use_eff_at_mod");
	}

	td->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT_TOT,"total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_E_PASS,"e^{-} ineff");
	hevtsum->GetXaxis()->SetBinLabel(EVT_P_PASS,"p ineff");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIP_PASS,"#pi^{+} ineff");
	hevtsum->GetXaxis()->SetBinLabel(EVT_E_AND_P_AND_PIP_PASS,"e+p+#pi^{+} ineff");
	hevtsum->SetMinimum(0.);

	//! cut monitor histograms
	//! _hthetavp[3][6][2] for e,p,pip:sector:prec,pstc
	_eff_prtcl_names=new TString[3];
	_eff_prtcl_names[0]="e";
	_eff_prtcl_names[1]="p";
	_eff_prtcl_names[2]="pip";
	_hthetavp=new TH2D***[3]; //!e,p,pip
	for (int i=0;i<3;i++){
		_hthetavp[i]=new TH2D**[6]; //!sector
		for (int j=0;j<6;j++){
			_hthetavp[i][j]=new TH2D*[2];//prec,pstc
			for (int k=0;k<2;k++){
				TString name_cut_lvl;
				if      (k==0) name_cut_lvl="prec";
				else if (k==1) name_cut_lvl="pstc";

				TString name=TString::Format("thetavp_%s_s%d_%s",
					                         _eff_prtcl_names[i].Data(),j+1,name_cut_lvl.Data());
				if      (_eff_prtcl_names[i]=="e")   _hthetavp[i][j][k]=new TH2D(name,name,100,0,5,100,0,60);
				else if (_eff_prtcl_names[i]=="p")   _hthetavp[i][j][k]=new TH2D(name,name,100,0,4,100,0,60);
				else if (_eff_prtcl_names[i]=="pip") _hthetavp[i][j][k]=new TH2D(name,name,100,0,3,100,0,120);	
			}
		}
	}

	
}

ProcEff::~ProcEff()
{
	delete[] _eff_prtcl_names;
	delete[] _hthetavp;
}

void ProcEff::handle()
{
	//Info("ProcEff::handle()", "");
	pass = kFALSE;
	hevtsum->Fill(EVT_TOT);

	updateEff();

	bool in_eff_rgn=kFALSE;
	if (dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0) { //! t1 and t2 pim does not matter; dAna->pidnew.h10IdxPim>0) {
		in_eff_rgn=pass_eff();
		if (in_eff_rgn) {
			hevtsum->Fill(EVT_E_AND_P_AND_PIP_PASS);
		}
	}else if (dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip==0 && dAna->pidnew.h10IdxPim>0){ //! t3
		in_eff_rgn=kTRUE;
	}else if (dAna->pidnew.h10IdxP==0 && dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0){//! t4
		in_eff_rgn=kTRUE;
	}
	
	if (in_eff_rgn)
	{
		pass = kTRUE;
		EpProcessor::handle();
	}
}

/* Currently implemented for top2'*/
bool ProcEff::pass_eff(){
	//Info("h10looper_2pi::pass_eff()","");
	bool ret=kFALSE;
	bool e_ineff=kFALSE;
	bool p_ineff=kFALSE;
	bool pip_ineff=kFALSE;

	DataEff* deff = &dAna->eff;

	int sector_e =deff->sector_e;
	float p_e    =deff->p_e;
	float theta_e=deff->theta_e;

	int sector_p =deff->sector_p;
	float p_p    =deff->p_p;
	float theta_p=deff->theta_p;

	int sector_pip =deff->sector_pip;
	float p_pip    =deff->p_pip;
	float theta_pip=deff->theta_pip;

	//Info("h10looper_2pi::pass_eff()","sectors:%d, %d, %d",_sector_e,_sector_p,_sector_pip);

	//!prec
	_hthetavp[0][sector_e-1][0]->Fill(p_e,theta_e);
	_hthetavp[1][sector_p-1][0]->Fill(p_p,theta_p);
	_hthetavp[2][sector_pip-1][0]->Fill(p_pip,theta_pip);

	//! check if e- ineff
	if (pass_theta_vs_p("e")){
		hevtsum->Fill(EVT_E_PASS);
		e_ineff=kTRUE;
		_hthetavp[0][sector_e-1][1]->Fill(p_e,theta_e);
	}
	//! check if p ineff
	if (pass_theta_vs_p("p")){
		hevtsum->Fill(EVT_P_PASS);
		p_ineff=kTRUE;
		_hthetavp[1][sector_p-1][1]->Fill(p_p,theta_p);
	}
	//! check if pip ineff
	if (pass_theta_vs_p("pip")){
		hevtsum->Fill(EVT_PIP_PASS);
		pip_ineff=kTRUE;
		_hthetavp[2][sector_pip-1][1]->Fill(p_pip,theta_pip);
	}

	ret=e_ineff && p_ineff && pip_ineff;
	return ret;
}

bool ProcEff::pass_theta_vs_p(TString prtcl_name){
	bool ret=kFALSE;
	DataEff* deff = &dAna->eff;
	
	int sector_e =deff->sector_e;
	float p_e    =deff->p_e;
	float theta_e=deff->theta_e;

	int sector_p =deff->sector_p;
	float p_p    =deff->p_p;
	float theta_p=deff->theta_p;

	int sector_pip =deff->sector_pip;
	float p_pip    =deff->p_pip;
	float theta_pip=deff->theta_pip;
	

	if (dH10->expt=="e1f"){
		ret=kTRUE;
	}else if (dH10->expt=="e16"){
		int at_mod=0;
		if (_use_eff_at_mod) at_mod=1;

		if (prtcl_name=="e") {
			ret=theta_vs_p_e16_el(sector_e,theta_e,p_e,at_mod);
		}else if (prtcl_name=="p") {
			ret=theta_vs_p_e16_pr(sector_p,theta_p,p_p,at_mod);
		}else if (prtcl_name=="pip") {
			ret=theta_vs_p_e16_pip(sector_pip,theta_pip,p_pip,at_mod);
		}
	}
	return ret;
}

void ProcEff::updateEff() {
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
	dAna->eff.sector_e=get_sector(phi_e);
	//! ******
	dAna->eff.theta_e=lvE.Theta()*RadToDeg();
	dAna->eff.p_e=lvE.P();

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
		dAna->eff.sector_p=get_sector(phi_p);
		//! ******
		dAna->eff.theta_p=lvP.Theta()*RadToDeg();
		dAna->eff.p_p=lvP.P();
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
		dAna->eff.sector_pip=get_sector(phi_pip);
		//! ******
		dAna->eff.theta_pip=lvPip.Theta()*RadToDeg();
		dAna->eff.p_pip=lvPip.P();
	}
}

/*
phi has to be [-30,330] for the following function to work.
*/
int ProcEff::get_sector(float phi){
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
		Info("ProcEff::get_sector","sector not found for phi=%.2f. Returning sector=-1",phi);
		//Info("ProcEff::get_sector","file=%s\n",fChain->GetCurrentFile()->GetName());
	}
	return sector;
}

#endif // PROCEFF_H
