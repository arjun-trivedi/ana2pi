#ifndef PROCDELAST_H
#define PROCDELAST_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "cuts.h"
#include "particle_constants.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

using namespace TMath;
using namespace ParticleConstants;

class ProcDelast : public EpProcessor {

public:
	ProcDelast(TDirectory *td,DataH10* dataH10,DataAna* dataAna,
			   bool monitor=kFALSE, bool monitorOnly=kFALSE,
		       bool procT=kFALSE, bool procR=kFALSE);
	~ProcDelast();
	void handle();
	void write();
protected:
	bool _procT,_procR;
	TDirectory* _dirmon;
	TDirectory* _dircut;
	//! For Thrown 
	TObjArray* _hists_ekin_ST;
	TTree* _t_ST;
	TObjArray* _yields_ST;
	//! For Reconstructed
	TObjArray* _hists_ekin[NPROCMODES];
	TObjArray* _hists_cuts[NPROCMODES];
	TTree* _t[NPROCMODES];
	TObjArray* _yields[NPROCMODES];
	bool _det_e,_det_p;
	int _h10idxE,_h10idxP;
	//! Used by ST,ER and SR
	TLorentzVector _lvE0;
	TLorentzVector _lvP0;
	TLorentzVector _lvQ;
	TLorentzVector _lvW;
	TLorentzVector _lvE;
	TLorentzVector _lvP;
	TLorentzVector _lvMMep;
	
	static const Int_t NUM_EVTCUTS = 2;
	enum { EVT_NULL, EVT,EVT_PASS};//, EVT_GPART, EVT_NE, EVT_NP};

	void AddBranches(TTree* t, bool ismc=kFALSE);

	void McKin();
	void UpdateDelast(bool ismc  = kFALSE);
	void updateEkin(Bool_t useMc = kFALSE);
	
	void ResetLvs();
	
	void setLabFrmLvs();
	bool passEvent();

	Float_t getPhi(TLorentzVector lv);
	Float_t invTan(Float_t y, Float_t x);


	//** Following 2 functions used to apply efid to ST-electrons
	bool inFid(); 
	float GetPhi(TLorentzVector lv);
	int GetSector(TLorentzVector lv); 
};

ProcDelast::ProcDelast(TDirectory *td,DataH10* dataH10,DataAna* dataAna,
					   bool monitor,bool monitorOnly,
	                   bool procT, bool procR
	                   ):EpProcessor(td, dataH10, dataAna,monitor,monitorOnly) {
	_procT=procT;
	_procR=procR;
	_dirmon=NULL;
	_dircut=NULL;
	_lvE0 = dH10->lvE0;
	_lvP0 = dH10->lvP0;

	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	//hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"nevts");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PASS,"passed");
}

ProcDelast::~ProcDelast() {
	//delete hevtsum;
	delete _dirmon,_dircut;
	
}

void ProcDelast::handle() {
	//Info("In ProcDelast::handle()");
	pass = kFALSE;
	hevtsum->Fill(EVT);

	ResetLvs();

	//! Only if _procT
	if (_procT){
		McKin();
		updateEkin(kTRUE);
		UpdateDelast(kTRUE);

		if (mon || mononly){//! If mon or mononly mode
			if (_t_ST==NULL){
				_t_ST = new TTree("t_ST","TTree containing data for Elastic events");
				AddBranches(_t_ST,kTRUE);
				_hists_ekin_ST=dAna->makeHistsEkin();
				_yields_ST=dAna->makeYieldsElastic();
			}
			_t_ST->Fill();
			dAna->fillHistsEkin(_hists_ekin_ST,kTRUE);
			dAna->fillYieldsElastic(_yields_ST,kTRUE);
		}else{//! If regular mode
			if (_hists_ekin_ST==NULL){
				_hists_ekin_ST=dAna->makeHistsEkin();
				_yields_ST=dAna->makeYieldsElastic();
			}
			dAna->fillHistsEkin(_hists_ekin_ST,kTRUE);
			dAna->fillYieldsElastic(_yields_ST,kTRUE);
		}
		pass=kTRUE;
		hevtsum->Fill(EVT_PASS);
		EpProcessor::handle(); 
		return;
	}

	//! Only if _procR
	if (_procR){
		setLabFrmLvs();
		updateEkin();
		UpdateDelast();

		if (mon || mononly){//! If mon or mononly mode
			if (_t[MONMODE]==NULL){
				Info("In ProcDelast::handle()","_t[MONMODE]!=NULL\n");
				_dirmon = dirout->mkdir(TString::Format("monitor"));
				_dirmon->cd();
				_t[MONMODE] = new TTree("t","TTree containing data for Elastic events");
				AddBranches(_t[MONMODE]);
				_hists_ekin[MONMODE]=dAna->makeHistsEkin();
				_hists_cuts[MONMODE]=dAna->makeHistsMMElastic();
				_yields[MONMODE]=dAna->makeYieldsElastic();
			}
			_t[MONMODE]->Fill();
			dAna->fillHistsEkin(_hists_ekin[MONMODE]);
			dAna->fillHistsMMElastic(_hists_cuts[MONMODE]);
			dAna->fillYieldsElastic(_yields[MONMODE]);
			if (mononly){
				pass=kTRUE;
				hevtsum->Fill(EVT_PASS);
				EpProcessor::handle();
				return;
			}
		}
		
		pass=passEvent(); //! Elastic event selection

		if (pass) {
			if (mon){//! If mon mode
				if (_t[CUTMODE]==NULL){
					Info("In ProcDelast::handle()","_t[CUTMODE]!=NULL\n");
					_dircut = dirout->mkdir(TString::Format("cut"));
					_dircut->cd();
					_t[CUTMODE] = new TTree("t","TTree containing data for Elastic events");
					AddBranches(_t[CUTMODE]);
					_hists_ekin[CUTMODE]=dAna->makeHistsEkin();
					_hists_cuts[CUTMODE]=dAna->makeHistsMMElastic();
					_yields[CUTMODE]=dAna->makeYieldsElastic();
				}
				_t[CUTMODE]->Fill();
				dAna->fillHistsEkin(_hists_ekin[CUTMODE]);
				dAna->fillHistsMMElastic(_hists_cuts[CUTMODE]);
				dAna->fillYieldsElastic(_yields[CUTMODE]);
			}else{//! If regular mode
				if (_hists_ekin[CUTMODE]==NULL){
					dirout->cd();
					_hists_ekin[CUTMODE]=dAna->makeHistsEkin();
					_hists_cuts[CUTMODE]=dAna->makeHistsMMElastic();
					_yields[CUTMODE]=dAna->makeYieldsElastic();
				}
				dAna->fillHistsEkin(_hists_ekin[CUTMODE]);
				dAna->fillHistsMMElastic(_hists_cuts[CUTMODE]);
				dAna->fillYieldsElastic(_yields[CUTMODE]);
			}
			hevtsum->Fill(EVT_PASS);
			EpProcessor::handle();
			return;
		}
	}
}

void ProcDelast::ResetLvs(){
	_lvQ.SetXYZT(0,0,0,0);
	_lvW.SetXYZT(0,0,0,0);
	_lvE.SetXYZT(0,0,0,0);
	_lvP.SetXYZT(0,0,0,0);
	_lvMMep.SetXYZT(0,0,0,0);
}

void ProcDelast::AddBranches(TTree* t, bool ismc/*=kFALSE*/){
	DataElastic *de = &(dAna->dElast);
	if (ismc) de = &(dAna->dElast_ST);

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
	t->Branch("nphe",&de->nphe);
}

void ProcDelast::setLabFrmLvs(){
	DataElastic *de = &(dAna->dElast);
	_det_e=kFALSE;
	_det_p=kFALSE;
	_h10idxE=-1;
	_h10idxP=-1;
		
	for (int igpart=0;igpart<dH10->gpart;igpart++)	{
		if (igpart==0 && dH10->id[igpart]==ELECTRON){
			_det_e=kTRUE;
			_h10idxE=igpart;	
			Double_t mom = dH10->p[igpart];
	       	Double_t px = mom*dH10->cx[igpart];
	       	Double_t py = mom*dH10->cy[igpart];
	       	Double_t pz = mom*dH10->cz[igpart];
	       	Double_t energy = Sqrt(mom*mom+MASS_E*MASS_E);
	       	_lvE.SetPxPyPzE(px,py,pz,energy);
		} 
		if (dH10->id[igpart]==PROTON){
			_det_p=kTRUE;
			_h10idxP=igpart;
			Double_t mom = dH10->p[igpart];
	       	Double_t px = mom*dH10->cx[igpart];
	       	Double_t py = mom*dH10->cy[igpart];
	       	Double_t pz = mom*dH10->cz[igpart];
	       	Double_t energy = Sqrt(mom*mom+MASS_P*MASS_P);
			_lvP.SetPxPyPzE(px,py,pz,energy);
		}   
	}
	//! Q2,W
	_lvQ = _lvE0-_lvE;
	_lvW = _lvQ+_lvP0;
	//! MMep
	_lvMMep=_lvQ+_lvP0-_lvP;

	return;
}

bool ProcDelast::passEvent(){
	bool ret=kFALSE;
	DataElastic *de = &(dAna->dElast);
	if (dH10->gpart<=4 && 
		_det_e==kTRUE && _det_p==kTRUE && 
		dH10->dc[_h10idxP]>0 && dH10->sc[_h10idxP]>0 &&
		de->MMep<0.1 && de->W<1.1)
		{
			ret=kTRUE;
		}
	return ret;
}

void ProcDelast::UpdateDelast(bool ismc /* = kFALSE */){
	DataElastic *de = &(dAna->dElast);
	if (ismc) de = &(dAna->dElast_ST);

	//gpart
	if (!ismc)	de->gpart=dH10->gpart;
	else 		de->gpart=dH10->nprt;
	//Q2,W
	de->Q2=-1*(_lvQ.Mag2());
	de->W=_lvW.Mag();
	//! MMs for Event selection
	de->MMep= _lvMMep.Mag();
	de->MM2ep=_lvMMep.Mag2();
	//! Reconstructed Kinematics 
	//! for e',p',p,pip,pim at e' vertex
	de->p_e=_lvE.P();
	de->p_p=_lvP.P();
	de->theta_e=_lvE.Theta()*RadToDeg();
	de->theta_p=_lvP.Theta()*RadToDeg();
	de->phi_e=getPhi(_lvE);//_lvE.Phi()*RadToDeg();
	de->phi_p=getPhi(_lvP);//_lvP.Phi()*RadToDeg();
	//! Reconstructed e' Vertex
	if (!ismc){
		de->vx_e=dH10->vx[_h10idxE];
		de->vx_p=dH10->vx[_h10idxP];
		de->vy_e=dH10->vy[_h10idxE];
		de->vy_p=dH10->vy[_h10idxP];
		de->vz_e=dH10->vz[_h10idxE];
		de->vz_p=dH10->vz[_h10idxP];
	}else{
		for (Int_t inprt=0; inprt<dH10->nprt;inprt++){
			if (dH10->pidpart[inprt]==3){
				de->vx_e=dH10->xpart[inprt];
				de->vy_e=dH10->ypart[inprt];
				de->vz_e=dH10->zpart[inprt];
			}
			if (dH10->pidpart[inprt]==14){
				de->vx_p=dH10->xpart[inprt];
				de->vy_p=dH10->ypart[inprt];
				de->vz_p=dH10->zpart[inprt];
			}
		}
	}
	//! EID
	if (!ismc){
		de->nphe=dH10->nphe[dH10->cc[_h10idxE]-1];
	}

}

void ProcDelast::McKin(){
	for (Int_t inprt=0; inprt<dH10->nprt;inprt++)	{
		if (dH10->pidpart[inprt]==3){
			Float_t px=dH10->pxpart[inprt];
			Float_t py=dH10->pypart[inprt];
			Float_t pz=dH10->pzpart[inprt];
			Float_t e=dH10->epart[inprt];
			_lvE.SetPxPyPzE(px,py,pz,e);
		}
		if (dH10->pidpart[inprt]==14){
			Float_t px=dH10->pxpart[inprt];
			Float_t py=dH10->pypart[inprt];
			Float_t pz=dH10->pzpart[inprt];
			Float_t e=dH10->epart[inprt];
			_lvP.SetPxPyPzE(px,py,pz,e);
		}
	}
	
	//! Q2,W
	_lvQ = _lvE0-_lvE;
	_lvW = _lvQ+_lvP0;
	//! MMep
	_lvMMep=_lvQ+_lvP0-_lvP;
}

bool ProcDelast::inFid() {
	TLorentzVector lvE(0,0,0,0);
	for (Int_t inprt=0; inprt<dH10->nprt;inprt++)	{
		if (dH10->pidpart[inprt]==3){
			Float_t px=dH10->pxpart[inprt];
			Float_t py=dH10->pypart[inprt];
			Float_t pz=dH10->pzpart[inprt];
			Float_t e=dH10->epart[inprt];
			lvE.SetPxPyPzE(px,py,pz,e);
		}
	}
	Int_t id = 11;
	Double_t p = lvE.P();
	Float_t theta = lvE.Theta()*RadToDeg();
	Float_t phi = lvE.Phi()*RadToDeg(); //phitmp = [-180,180]
	Int_t sector=GetSector(lvE);
	Int_t paddle=0;
	bool inFid = Cuts::Fiducial(id,p,theta,phi,sector,paddle);
	return inFid;
}

int ProcDelast::GetSector(TLorentzVector lv){
	float phi = GetPhi(lv);

	int sector = -9999;
	if( (330 < phi && phi <= 360) || (0 <= phi && phi <= 30) ){
		sector = 1;
	}
	if(30 < phi && phi <= 90){
		sector = 2;
	}
	if(90 < phi && phi <= 150){
		sector = 3;
	}
	if(150 < phi && phi <= 210){
		sector = 4;
	}
	if(210 < phi && phi <= 270){
		sector = 5;
	}
	if(270 < phi && phi <= 330){
		sector = 6;
	}
	return sector;
}

float ProcDelast::GetPhi(TLorentzVector lv){
	double phi = -9999;
	/* TLorentzVector Phi() returns angle between -pi and pi using ATan2(y,x)
	   The following transforms [-pi,pi]-->[0,2pi] */
	if(lv.Phi()*(180/TMath::Pi()) < 0){
		phi = lv.Phi()*(180/TMath::Pi()) + 360;
	} else {
		phi = lv.Phi()*(180/TMath::Pi());
	}
	return phi;
}

void ProcDelast::updateEkin(Bool_t useMc /*= kFALSE*/) {
	const TLorentzVector _4vE0 = dH10->lvE0;
	const TLorentzVector _4vP0 = dH10->lvP0;
	TLorentzVector _4vE1;
		
	DataEkin *ekin = &dAna->eKin;
	if (useMc) ekin = &dAna->eKin_mc;
	
	Float_t mom;
	Float_t px;
	Float_t py;
	Float_t pz;
	if (!useMc){
		mom = dH10->p[0];
		px = mom*dH10->cx[0];
		py = mom*dH10->cy[0];
		pz = mom*dH10->cz[0];
		_4vE1.SetPxPyPzE(px,py,pz,Sqrt(mom*mom+MASS_E*MASS_E));
	} else{
		for (Int_t inprt=0; inprt<dH10->nprt;inprt++) {
			if (dH10->pidpart[inprt]==3){
				Float_t px=dH10->pxpart[inprt];
				Float_t py=dH10->pypart[inprt];
				Float_t pz=dH10->pzpart[inprt];
				Float_t e=dH10->epart[inprt];
				_4vE1.SetPxPyPzE(px,py,pz,e);
			}
		}
	}
	
	TLorentzVector _4vQ = _4vE0-_4vE1;
	if (!useMc) {ekin->sector = dH10->sc_sect[dH10->sc[0]-1];}
	ekin->W = (_4vQ+_4vP0).Mag();
	ekin->Q2 = -1*_4vQ.Mag2();
	ekin->nu = _4vQ.E();
	ekin->xb = ekin->Q2/(2*MASS_P*ekin->nu);
	ekin->E1 = _4vE1.E();
	ekin->theta1 = _4vE1.Theta()*RadToDeg();
	Double_t phitmp = _4vE1.Phi()*RadToDeg(); //phitmp = [-180,180]
	ekin->phi1 = phitmp<-30 ? phitmp+360 : phitmp;
	ekin->theta = _4vQ.Theta()*RadToDeg();
	phitmp = _4vQ.Phi()*RadToDeg(); //phitmp = [-180,180]
	ekin->phi = phitmp<-30 ? phitmp+360 : phitmp;
}

Float_t ProcDelast::getPhi(TLorentzVector lv) {
	Float_t retVal = 0;
	retVal = RadToDeg()*invTan(lv.Py(), lv.Px());
	return retVal;
}

Float_t ProcDelast::invTan(Float_t y, Float_t x){
	Float_t retVal = 0; 
	if (x > 0 && y > 0)  retVal = ATan(y/x);          //1st Quad.
	if (x < 0 && y > 0)  retVal = ATan(y/x) + Pi();   //2nd Quad
	if (x < 0 && y < 0)  retVal = ATan(y/x) + Pi();   //3rd Quad
	if (x > 0 && y < 0)  retVal = ATan(y/x) + 2*Pi(); //4th Quad 
	if (x == 0 && y > 0) retVal = Pi()/2;
	if (x == 0 && y < 0) retVal = 3*Pi()/2; 
	return retVal;  
}


void ProcDelast::write(){

}

#endif // PROCDELAST_H
