#ifndef PROCEID_H
#define PROCEID_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "eid.h"
#include "particle_constants.h"
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TMath.h>

using namespace TMath;
using namespace ParticleConstants;
using namespace AnalysisConstants;

class ProcEid : public EpProcessor
{
public:
	ProcEid(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		    Bool_t mon=kFALSE,Bool_t monOnly=kFALSE);
	~ProcEid();
	
	void handle();
	//void write();
		
protected:
	Eid* _eidTool;
				
	static const Int_t NUM_EVTCUTS = 12;
	enum { EVT_NULL, EVT_TRIG, EVT_GPART1, EVT_STAT1, EVT_Q1,
	       EVT_SC1, EVT_DC1, EVT_EC1, EVT_CC1,
	       EVT_DCSTAT1, EVT_ECLOW1, EVT_SF, EVT_BOS11
	     };
	
	Bool_t goodE(DataH10* dH10);
	Bool_t goodE_bos(DataH10* dH10);
	void updateEid(DataH10* dH10);
	void updateEkin(DataH10* dH10, Bool_t useMc = kFALSE);
	Float_t getCCtheta(Float_t x_sc, Float_t y_sc, Float_t z_sc, Float_t cx_sc, Float_t cy_sc, Float_t cz_sc);
};

ProcEid::ProcEid(TDirectory *td, DataH10* dataH10, DataAna* dataAna, 
                 Bool_t mon/* = kFALSE*/,Bool_t monOnly /*= kFALSE*/)
                 :EpProcessor(td, dataH10, dataAna, mon, monOnly)
{
	if      (dH10->expt=="e1f" && dH10->dtyp=="sim") _eidTool = new Eid("/home/trivedia/CLAS/workspace/ana2pi/eid/eid.mc.out");
	else if (dH10->expt=="e1f" && dH10->dtyp=="exp") _eidTool = new Eid("/home/trivedia/CLAS/workspace/ana2pi/eid/eid.exp.out");
	else  Info("ProcEid::ProcEid()", "_eidTool not initialized");//for e1-6

    if      (dH10->expt=="e1f" && _eidTool->eidParFileFound) {
    	Info("ProcEid::ProcEid()", "is_h10e1f=true && eidParFileFound=true. Will use goodE()"); 
    }else if (dH10->expt=="e1f" && !_eidTool->eidParFileFound) {
    	Info("ProcEid::ProcEid()", "is_h10e1f=true && eidParFileFound=false. Will use goodE_bos()");
    }else if (dH10->expt=="e16") {
    	Info("ProcEid::ProcEid()", "is_h10e16=true. Will use goodE_bos()");; //pars for e1-6 not yet obtained
    }
	
	td->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	//hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT_TRIG,"Trigger");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPART1,"gpart>1");
	hevtsum->GetXaxis()->SetBinLabel(EVT_STAT1,"stat>0");
	hevtsum->GetXaxis()->SetBinLabel(EVT_Q1,"q=-1");
	hevtsum->GetXaxis()->SetBinLabel(EVT_SC1,"SC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_DC1,"DC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_EC1,"EC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ECLOW1,"EC Threshold");
	hevtsum->GetXaxis()->SetBinLabel(EVT_CC1,"CC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_DCSTAT1,"Time-based");
	hevtsum->GetXaxis()->SetBinLabel(EVT_SF,"SF");
	hevtsum->GetXaxis()->SetBinLabel(EVT_BOS11,"EVNT.id=11");
}

ProcEid::~ProcEid(){
	delete _eidTool;
}
	
void ProcEid::handle() {
	//Info("ProcEid::handle()", "");
	pass = kFALSE;
	
	hevtsum->Fill(EVT_TRIG);
	
	if (mMon||mMonOnly)
	{
		if (dAna->top==0 && hists[iMODE_MON][iEVTINC][iSECTOR0]==NULL) { //i.e. inclusive event
			TDirectory* dirmon = dirout->mkdir(TString::Format("mon"));
			dAna->makeHistsEid(hists[iMODE_MON][iEVTINC], dirmon);
			dAna->makeHistsEkin(histsEkin[iMODE_MON][iEVTINC], dirmon);
		}else if(dAna->top!=0 && hists[iMODE_MON][iTOP1][iSECTOR0]==NULL){ //i.e. 2pi event
			for(Int_t iTop=iTOP1;iTop<nTOP;iTop++){
				TDirectory* dirmon = dirout->mkdir(TString::Format("mon%d",iTop));
				dAna->makeHistsEid(hists[iMODE_MON][iTop], dirmon);
				dAna->makeHistsEkin(histsEkin[iMODE_MON][iTop], dirmon);
			}
		}
		
		updateEid(dH10);
		updateEkin(dH10);
	
		//if ( (dAna->eKin.W < 1.5) ) return;
		
		if (dAna->top == 0) { //i.e inclusive event
			dAna->fillHistsEid(hists[iMODE_MON][iEVTINC]);
			dAna->fillHistsEkin(histsEkin[iMODE_MON][iEVTINC]);
		}else { //i.e 2pi event
			dAna->fillHistsEid(hists[iMODE_MON][dAna->top]);
			dAna->fillHistsEkin(histsEkin[iMODE_MON][dAna->top]);
		}
	}
	
		
	if (mMonOnly){
		pass = kTRUE;
		EpProcessor::handle();
		return;
	}
	
	Bool_t gE = kFALSE;

    if      (is_h10e1f && _eidTool->eidParFileFound)  gE =  goodE(dH10);
    else if (is_h10e1f && !_eidTool->eidParFileFound) gE =  goodE_bos(dH10);
    else if (is_h10e16)                               gE =  goodE_bos(dH10); //pars for e1-6 not yet obtained
    
	if (gE) {	
		
		if (mMon)
		{
			if (hists[iMODE_CUT][iEVTINC][iSECTOR0]==NULL) {
				TDirectory* dircut = dirout->mkdir(TString::Format("cut"));
				dAna->makeHistsEid(hists[iMODE_CUT][iEVTINC], dircut);
				dAna->makeHistsEkin(histsEkin[iMODE_CUT][iEVTINC], dircut);
			}
			dAna->fillHistsEid(hists[iMODE_CUT][iEVTINC]);
			dAna->fillHistsEkin(histsEkin[iMODE_CUT][iEVTINC]);
		}
		
		dH10->id[0] = ELECTRON;
		pass = kTRUE;
		EpProcessor::handle();
	}
}

Bool_t ProcEid::goodE(DataH10* dH10){
	Bool_t retval = kFALSE;
	
	if (dH10->dtyp!="sim") { //atrivedi 020313 till _eidTool is fixed to use eid.mc.out
		if (dH10->id[0]==ELECTRON) hevtsum->Fill(EVT_BOS11);
	}
	if (dH10->gpart>1) {
		hevtsum->Fill(EVT_GPART1);
		if (dH10->stat[0]>0) {
			hevtsum->Fill(EVT_STAT1);
			if (dH10->q[0]==-1) {
				hevtsum->Fill(EVT_Q1);
				if (dH10->sc[0]>0) {
					hevtsum->Fill(EVT_SC1);
					if (dH10->dc[0]>0) {
						hevtsum->Fill(EVT_DC1);
						if (dH10->ec[0]>0) {
							hevtsum->Fill(EVT_EC1);
							if (dH10->cc[0]>0 || dH10->dtyp=="sim") {
								hevtsum->Fill(EVT_CC1);
								if (dH10->dc_stat[dH10->dc[0]-1]>0) {
									hevtsum->Fill(EVT_DCSTAT1);
									if (_eidTool->PassThreshold(dH10->p[0])) {
										hevtsum->Fill(EVT_ECLOW1);
										Int_t sector = dH10->sc_sect[dH10->sc[0]-1];
										Float_t mom = dH10->p[0];
										Float_t sf = dH10->etot[dH10->ec[0]-1]/dH10->p[0];
										if (_eidTool->PassSF(sector,mom,sf)) {
											hevtsum->Fill(EVT_SF);
											if (!_eidTool->Pass(sector,mom,sf)) {
												hevtsum->Fill(NUM_EVTCUTS+1);
											}
											retval = kTRUE;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return retval;
}

Bool_t ProcEid::goodE_bos(DataH10* dH10){
	Bool_t retval = kFALSE;
	
	//if (dH10->gpart>1) {
		//hevtsum->Fill(EVT_GPART1);
		//if (dH10->stat[0]>0) {
			//hevtsum->Fill(EVT_STAT1);
			if (dH10->q[0]==-1) {
				hevtsum->Fill(EVT_Q1);
				if (dH10->sc[0]>0) {
					hevtsum->Fill(EVT_SC1);
					if (dH10->dc[0]>0) {
						hevtsum->Fill(EVT_DC1);
						if (dH10->ec[0]>0) {
							hevtsum->Fill(EVT_EC1);
							if (dH10->cc[0]>0 || dH10->dtyp=="sim") {
								hevtsum->Fill(EVT_CC1);
								//if (dH10->dc_stat[dH10->dc[0]-1]>0) {
									//hevtsum->Fill(EVT_DCSTAT1);
									if (dH10->id[0]==ELECTRON) {
										hevtsum->Fill(EVT_BOS11);
										retval = kTRUE;
									}
								//}
							}
						}
					}
				}
			}
		//}
	//}
	return retval;
}

void ProcEid::updateEid(DataH10* dH10){
		
	Float_t p = dH10->p[0];
	Float_t l_e = dH10->sc_r[dH10->sc[0]-1];
	Float_t t_e = dH10->sc_t[dH10->sc[0]-1];
	Float_t tOFF = t_e-(l_e/SOL);
	
	Float_t beta = ( l_e/(t_e - tOFF) )/SOL; // = 1
	Float_t betaStrE = TMath::Sqrt((p*p)/(MASS_E*MASS_E+p*p));
	Float_t dtE = l_e/(betaStrE*SOL) + tOFF - t_e; // = 0
	
	dAna->eid.sector = dH10->sc_sect[dH10->sc[0]-1];
	
	dAna->eid.beta = beta;
	dAna->eid.betaStrE = betaStrE;
	dAna->eid.dtE = dtE;
		
	dAna->eid.p = p;
	dAna->eid.sector  = dH10->sc_sect[dH10->sc[0]-1];
	//from EC
	dAna->eid.ec_ei   = dH10->ec_ei[dH10->ec[0]-1];
	dAna->eid.ec_eo   = dH10->ec_eo[dH10->ec[0]-1];
	dAna->eid.etot    = dH10->etot[dH10->ec[0]-1];
	//from CC
	dAna->eid.nphe    = dH10->nphe[dH10->cc[0]-1];
	dAna->eid.cc_segm = (dH10->cc_segm[dH10->cc[0]-1]%1000)/10;
	//calculate cc_theta
	Float_t dc_xsc  = dH10->dc_xsc[dH10->dc[0]-1];
	Float_t dc_ysc  = dH10->dc_ysc[dH10->dc[0]-1];
	Float_t dc_zsc  = dH10->dc_zsc[dH10->dc[0]-1];
	Float_t dc_cxsc = dH10->dc_cxsc[dH10->dc[0]-1];
	Float_t dc_cysc = dH10->dc_cysc[dH10->dc[0]-1];
	Float_t dc_czsc = dH10->dc_czsc[dH10->dc[0]-1];
	dAna->eid.cc_theta = getCCtheta(dc_xsc, dc_ysc, dc_zsc, dc_cxsc, dc_cysc, dc_czsc);
}

void ProcEid::updateEkin(DataH10* dH10, Bool_t useMc /*= kFALSE*/) {
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
		for (Int_t idx = 0; idx < dH10->mcnentr; idx++) {
			Int_t _id = dH10->mcid[idx];
			if (_id != ELECTRON) continue;
			mom = dH10->mcp[idx];
			Float_t _theta = dH10->mctheta[idx]*DegToRad();
			Float_t _phi = dH10->mcphi[idx]*DegToRad();
			px = mom*Cos(_phi)*Sin(_theta);
			py = mom*Sin(_phi)*Sin(_theta);
			pz = mom*Cos(_theta);
			_4vE1.SetPxPyPzE(px,py,pz,Sqrt(mom*mom+MASS_E*MASS_E));
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

Float_t ProcEid::getCCtheta(Float_t x_sc, Float_t y_sc, Float_t z_sc, Float_t cx_sc, Float_t cy_sc, Float_t cz_sc){
	//Define CC plane equation: Ax + By + Cz + D = 0
	Float_t A = -0.000785;
	Float_t B = 0;
	Float_t C = -0.00168;
	Float_t D = 1;
	Float_t magS = TMath::Sqrt(A*A + B*B + C*C);

	//length of line perpendicular to CC plane from hit in SC plane
	Float_t h = 0;

	//cosine of angle between h & t vector
	Float_t cAlpha = 0;

	//magnitude of vector t
	Float_t t = 0;

	h = ((A*x_sc + B*y_sc + C*z_sc) + D)/magS;

	cAlpha = (cx_sc*x_sc + cy_sc*y_sc + cz_sc*z_sc)/magS;

	t = h/cAlpha;

	Float_t x_cc = x_sc - t*cx_sc;
	Float_t y_cc = y_sc - t*cy_sc;
	Float_t z_cc = z_sc - t*cz_sc;

	return TMath::ACos(z_cc/TMath::Sqrt(x_cc*x_cc + y_cc*y_cc + z_cc*z_cc)) * (180/TMath::Pi());
}
#endif // PROCEID_H
