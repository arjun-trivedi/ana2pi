#ifndef PROCPIDELAST_H
#define PROCPIDELAST_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"

using namespace ParticleConstants;

/******************************************************
[03-08-15]
+ Note, ProcPidElast must follow ProcSkimQElast
+ Currently using BOS PID
+ Pass Events only if:
	+ dAna->skimq_elast.isEVT_1POS_EX AND proton is identified
*******************************************************/

class ProcPidElast : public EpProcessor
{

public:
	ProcPidElast(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		    Bool_t monitor=kFALSE,Bool_t monitorOnly=kFALSE);
	ProcPidElast(DataH10* dataH10,DataAna* dataAna);
	~ProcPidElast();

	void handle();
	//void write();
	void updatePid();
				
protected:
	static const Int_t NUM_EVTCUTS=3;
	enum { EVT_NULL, EVT, EVT_P_EX, EVT_OTHER};   
	     
	//void updatePid();
	Float_t getCCtheta(Float_t x_sc, Float_t y_sc, Float_t z_sc, Float_t cx_sc, Float_t cy_sc, Float_t cz_sc);
};

ProcPidElast::ProcPidElast(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
                 Bool_t monitor/* = kFALSE*/,Bool_t monitorOnly /*= kFALSE*/)
                 :EpProcessor(td, dataH10, dataAna, monitor, monitorOnly)
{
	td->cd();	
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_P_EX,"p");
	hevtsum->GetXaxis()->SetBinLabel(EVT_OTHER,"other");
}

ProcPidElast::ProcPidElast(DataH10* dataH10,DataAna* dataAna)
                 :EpProcessor(dataH10, dataAna)
{
	
}

ProcPidElast::~ProcPidElast()
{
	
}

void ProcPidElast::handle()
{
	//Info("ProcPidElast::handle()", "");
	pass = kFALSE;
	
	hevtsum->Fill(EVT);
	
	if (mon||mononly)
	{
		if (hists[MONMODE][EVTINC][SECTOR0]==NULL) { //i.e. inclusive event
			TDirectory* dirmon = dirout->mkdir(TString::Format("monitor"));
			dAna->makeHistsPidElast(hists[MONMODE][EVTINC], dirmon);
		}	
	
		updatePid();
		if (dAna->d2pi.top == 0) { //i.e inclusive event
			dAna->fillHistsPidElast(hists[MONMODE][EVTINC]);
		}
	}
	
		
	if (mononly){
		pass = kTRUE;
		EpProcessor::handle();
		return;
	}	
	
	//pid cut
	dAna->pid_elast.h10IdxP=0; //atrivedi 041713: First reset indices set in monitoring mode; better soln?
	for (Int_t i = 1; i < dH10->gpart; i++) {
		if(dH10->q[i]==1){
			if (dH10->dc[i]>0) {
				if (dH10->sc[i]>0) {
					if (dH10->id[i]==PROTON){
						dH10->id[i]=PROTON;
						dAna->pid_elast.h10IdxP=i;
					}
				}
			}
		}
	}
	if(dAna->skimq_elast.isEVT_1POS_EX){
	 	if (dAna->pid_elast.h10IdxP>0) {
			hevtsum->Fill(EVT_P_EX);
			pass = kTRUE;
		}else{
			hevtsum->Fill(EVT_OTHER);
		}
	}
		
	if (pass) {
		if (mon)
		{
			if (hists[CUTMODE][EVTINC][SECTOR0]==NULL) {
				TDirectory* dircut = dirout->mkdir(TString::Format("cut"));
				dAna->makeHistsPidElast(hists[CUTMODE][EVTINC], dircut);
			}
			dAna->fillHistsPidElast(hists[CUTMODE][EVTINC]);
		}
		EpProcessor::handle();
	}
}

void ProcPidElast::updatePid()
{
	Float_t l_e = dH10->sc_r[dH10->sc[0]-1];
	Float_t t_e = dH10->sc_t[dH10->sc[0]-1];
	Float_t tOFF = t_e-(l_e/SOL);
	
	for (Int_t i = 1; i < dH10->gpart; i++) {
		Float_t p = dH10->p[i];
		Float_t l = dH10->sc_r[dH10->sc[i]-1];
		Float_t t = dH10->sc_t[dH10->sc[i]-1];
		Float_t beta = ( l/(t - tOFF) )/SOL; 
					
		if(dH10->id[i] == PROTON){
			Float_t betaStrP = TMath::Sqrt((p*p)/(MASS_P*MASS_P+p*p));
			Float_t dtP = l/(betaStrP*SOL) + tOFF - t;
			
			dAna->pid_elast.h10IdxP = i;
			dAna->pid_elast.sectorP = dH10->sc_sect[dH10->sc[i]-1];
			dAna->pid_elast.betaP = beta;
			dAna->pid_elast.pP = p;
			dAna->pid_elast.betaStrP = betaStrP;
			dAna->pid_elast.dtP = dtP;
			//from EC
			dAna->pid_elast.P_ec_ei   = dH10->ec_ei[dH10->ec[i]-1];
			dAna->pid_elast.P_ec_eo   = dH10->ec_eo[dH10->ec[i]-1];
			dAna->pid_elast.P_etot    = dH10->etot[dH10->ec[i]-1];
			//from CC
			dAna->pid_elast.P_nphe    = dH10->nphe[dH10->cc[i]-1];
			dAna->pid_elast.P_cc_segm = (dH10->cc_segm[dH10->cc[i]-1]%1000)/10;
			//calculate cc_theta
			Float_t dc_xsc  = dH10->dc_xsc[dH10->dc[i]-1];
			Float_t dc_ysc  = dH10->dc_ysc[dH10->dc[i]-1];
			Float_t dc_zsc  = dH10->dc_zsc[dH10->dc[i]-1];
			Float_t dc_cxsc = dH10->dc_cxsc[dH10->dc[i]-1];
			Float_t dc_cysc = dH10->dc_cysc[dH10->dc[i]-1];
			Float_t dc_czsc = dH10->dc_czsc[dH10->dc[i]-1];
			dAna->pid_elast.P_cc_theta = getCCtheta(dc_xsc, dc_ysc, dc_zsc, dc_cxsc, dc_cysc, dc_czsc);
		}
	}
}

Float_t ProcPidElast::getCCtheta(Float_t x_sc, Float_t y_sc, Float_t z_sc, Float_t cx_sc, Float_t cy_sc, Float_t cz_sc){
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

#endif // PROCPIDELAST_H
