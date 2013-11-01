#include "proc_tops.h"
#include "particle_constants.h"
#include "cuts.h"
#include <TMath.h>

using namespace TMath;
using namespace ParticleConstants;

ProcTops::ProcTops(TDirectory *td) : EpProcessor(td) {
	//histsMM = NULL; //for techincal reasons, hists was used
	histseKin = NULL;
	yields_top0 = NULL;
	yields_top1 = NULL;
	yields_top2 = NULL;
	yields_top3 = NULL;
	yields_top4 = NULL;
	histsMM_top0 = NULL;
	histsMM_top1 = NULL;
	histsMM_top2 = NULL;
	histsMM_top3 = NULL;
	histsMM_top4 = NULL;
	histseKin_top0 = NULL;
	histseKin_top1 = NULL;
	histseKin_top2 = NULL;
	histseKin_top3 = NULL;
	histseKin_top4 = NULL;  
	yields_mc = NULL;
	histsMM_mc = NULL;
	histseKin_mc = NULL; 
	dirout->cd();
	hevtsum = new TH1F("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPART0,"gpart > 0");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPARTEQ1,"gpart = 1");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPARTEQ2,"gpart = 2");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPARTEQ3,"gpart = 3");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPARTEQ4,"gpart = 4");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPART4,"gpart > 4");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GOODE,"good e");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GOODP,"good p");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GOODPIP,"good pi^{+}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GOODPIM,"good pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_TOP0,"Top0"); //TOPO = PPIPPIM||PPIP||PPIM||PIPPIM
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIPPIM,"Top1(p#pi^{+}#pi^{-})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIP,"Top2(p#pi^{+})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIM,"Top3(p#pi^{-})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIPPIM,"Top4(#pi^{+}#pi^{-})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_OTHER,"other");
	//hevtsum->SetDirectory(dirout);
}

ProcTops::~ProcTops() {
}

void ProcTops::handle(DataOmega* d) {
	//printf("In ProcTops::handle\n");
	hevtsum->Fill(EVT);
	pass = kFALSE;
	lvQ.SetXYZT(0,0,0,0);
	lvW.SetXYZT(0,0,0,0);
	lvE.SetXYZT(0,0,0,0);
	lvP.SetXYZT(0,0,0,0);
	lvPip.SetXYZT(0,0,0,0);
	lvPim.SetXYZT(0,0,0,0);
	lvMMtop1.SetXYZT(0,0,0,0);
	lvMMtop2.SetXYZT(0,0,0,0);
	lvMMtop3.SetXYZT(0,0,0,0);
	lvMMtop4.SetXYZT(0,0,0,0);
	if (hists==NULL) {
		dirout->cd();
		hists = d->makeHistsMM();
		histseKin = d->makeHistsEkin();
		dirout->mkdir("top0")->cd();
		yields_top0 = d->makeYields();
		histsMM_top0  = d->makeHistsMM();
		histseKin_top0 = d->makeHistsEkin();
		dirout->mkdir("top1")->cd();
		yields_top1 = d->makeYields();
		histsMM_top1  = d->makeHistsMM();
		histseKin_top1 = d->makeHistsEkin();
		dirout->mkdir("top2")->cd();
		yields_top2 = d->makeYields();
		histsMM_top2  = d->makeHistsMM();
		histseKin_top2 = d->makeHistsEkin();
		dirout->mkdir("top3")->cd();
		yields_top3 = d->makeYields();
		histsMM_top3  = d->makeHistsMM();
		histseKin_top3 = d->makeHistsEkin();
		dirout->mkdir("top4")->cd();
		yields_top4 = d->makeYields();
		histsMM_top4  = d->makeHistsMM();
		histseKin_top4 = d->makeHistsEkin();
		if (d->h10.is_sim) {
			dirout->cd();
			dirout->mkdir("mc")->cd();
			yields_mc = d->makeYields();
			histsMM_mc  = d->makeHistsMM();
			histseKin_mc = d->makeHistsEkin();
		}
	}
	if (d->h10.is_sim) {
		McKin(d);
		d->fillYields(yields_mc, kTRUE);
		d->fillHistsMM(histsMM_mc, kTRUE);
		d->fillHistsEkin(histseKin_mc, kTRUE);
	}
	if (d->h10.gpart>0) {hevtsum->Fill(EVT_GPART0);}
	if (d->h10.gpart==1) {hevtsum->Fill(EVT_GPARTEQ1);}
	if (d->h10.gpart==2) {hevtsum->Fill(EVT_GPARTEQ2);}
	if (d->h10.gpart==3) {hevtsum->Fill(EVT_GPARTEQ3);}
	if (d->h10.gpart==4) {hevtsum->Fill(EVT_GPARTEQ4);}
	if (d->h10.gpart>4) {hevtsum->Fill(EVT_GPART4);}
	for (Int_t i = 0; i < d->h10.gpart; i++) {
		Int_t id = d->h10.id[i];
		Double_t p = d->h10.p[i];
		Float_t theta = RadToDeg()*ACos(d->h10.cz[i]);
		Float_t phi = RadToDeg()*ATan2(d->h10.cy[i],d->h10.cx[i]);
		Int_t scidx = d->h10.sc[i]-1;
		Int_t sector = 1;
		Int_t paddle = 0;
		if (scidx >= 0) {
			sector = d->h10.sc_sect[scidx];
			paddle = d->h10.sc_pd[scidx];
		}
		Bool_t inFid = Cuts::Fiducial(id, p, theta, phi, sector, paddle);
		switch (id) {
		case ELECTRON:
			d->h10idxE = i;
			d->fid.fidE = inFid;
			//hevtsum->Fill(EVT_GOODE);
			break;
		case PROTON:
			if (d->h10.q[i] == 1 && d->h10.dc[i] > 0 && d->h10.sc[i] > 0) {
				d->h10idxP = i;
			    d->fid.fidP = inFid;
			    //hevtsum->Fill(EVT_GOODP);
			}
			break;
		case PIP:
			if (d->h10.q[i] == 1 && d->h10.dc[i] > 0 && d->h10.sc[i] > 0) {
				d->h10idxPip = i;
			    d->fid.fidPip = inFid;
			    //hevtsum->Fill(EVT_GOODPIP);
			}
			break;
		case PIM:
		    if (d->h10.q[i] == -1 && d->h10.dc[i] > 0 && d->h10.sc[i] > 0) {
				d->h10idxPim = i;
			    d->fid.fidPim = inFid;
			    //hevtsum->Fill(EVT_GOODPIM);
			}
			break;
		default:
			break;
		}
	}
	Bool_t gE = d->h10idxE>-1 && d->fid.fidE; //e1f
	if (gE) {hevtsum->Fill(EVT_GOODE);} 
	Bool_t gP = d->h10idxP>0;     //(no FidCut for p yet)
	if (gP) {hevtsum->Fill(EVT_GOODP);}
	Bool_t gPip = d->h10idxPip>0; //(no FidCut for p yet)
	if (gPip) {hevtsum->Fill(EVT_GOODPIP);}
	Bool_t gPim = d->h10idxPim>0; //(no FidCut for pim yet)
	if (gPim) {hevtsum->Fill(EVT_GOODPIM);}
	if (gE) { 
		if ( (gP && gPip && gPim) || (gP && gPip) || (gP && gPim) || (gPip && gPim) ) {
			/* *** Q2, W *** */
			Double_t mom = d->h10.p[d->h10idxE];
	        Double_t px = mom*d->h10.cx[d->h10idxE];
	        Double_t py = mom*d->h10.cy[d->h10idxE];
	        Double_t pz = mom*d->h10.cz[d->h10idxE];
	        Double_t energy = Sqrt(mom*mom+MASS_E*MASS_E);
	        lvE.SetPxPyPzE(px,py,pz,energy);
	        lvQ = _4vE0-lvE;
	        lvW = lvQ+_4vP0;
	        d->twoPi.Q2 = -1*(lvQ.Mag2());
	        d->twoPi.W = lvW.Mag();
	        
	        /* *** lvP, lvPip, lvPim *** */
	        set4vecHadrons(d);
	              
	        /* *** MM *** */         
			calcMM(d);
			Float_t mmppippim = Sqrt(d->twoPi.mm2ppippim);
			Float_t mmppip = Sqrt(d->twoPi.mm2ppip);
			Float_t mmppim = Sqrt(d->twoPi.mm2ppim);
			Float_t mmpippim = Sqrt(d->twoPi.mm2pippim);
			//Bool_t pi_mass = mmp > 0.7 && mmp < 0.863;
			
			Bool_t t1a = gP && gPip && gPim && d->h10.gpart == 4;
			Bool_t t1b = TMath::Abs(d->twoPi.mm2ppippim) < 0.0005;
			Bool_t t2a = gP && gPip && !gPim && d->h10.gpart == 3;
			Bool_t t2b = d->twoPi.mm2ppip>0 && d->twoPi.mm2ppip<0.04;
			Bool_t t3a = gP && gPim && !gPip && d->h10.gpart == 3;
			Bool_t t3b = d->twoPi.mm2ppim>0 && d->twoPi.mm2ppim<0.04;
			Bool_t t4a = gPip && gPim && !gP && d->h10.gpart == 3;
			Bool_t t4b = d->twoPi.mm2pippim>0.8 && d->twoPi.mm2pippim<1;
			
			Bool_t t1 = t1a && t1b;
			Bool_t t2 = t2a && t2b;
			Bool_t t3 = t3a && t3b;
			Bool_t t4 = t4a && t4b;
			
			if (t1a || t2a || t3a || t4a) {
				//calcVars(d);
				d->updateEkin();
				d->fillHistsMM(hists);
				d->fillHistsEkin(histseKin);
			}
			if ( t1 || t2 || t3 || t4) {
				pass = kTRUE;
				
				hevtsum->Fill(EVT_TOP0);
				d->topology = 0;
				calcVars(d);
				d->fillYields(yields_top0);
				d->fillHistsMM(histsMM_top0);
				d->fillHistsEkin(histseKin_top0);
				
				if (t1) {
					hevtsum->Fill(EVT_PPIPPIM);
					d->topology = 1;
					//calcVars(d);
					d->fillYields(yields_top1);
					d->fillHistsMM(histsMM_top1);
					d->fillHistsEkin(histseKin_top1);
				}
				if (t2) {
					hevtsum->Fill(EVT_PPIP);
					d->topology = 2;
					lvPim = lvMMtop2;
					//calcVars(d);
					d->fillYields(yields_top2);
					d->fillHistsMM(histsMM_top2);
					d->fillHistsEkin(histseKin_top2);
				}
				if (t3) {
					hevtsum->Fill(EVT_PPIM);
					d->topology = 3;
					lvPip = lvMMtop3;
					//calcVars(d);
					d->fillYields(yields_top3);
					d->fillHistsMM(histsMM_top3);
					d->fillHistsEkin(histseKin_top3);
				}
				if (t4) {
					hevtsum->Fill(EVT_PIPPIM);
					d->topology = 4;
					lvP = lvMMtop4;
					//calcVars(d);
					d->fillYields(yields_top4);
					d->fillHistsMM(histsMM_top4);
					d->fillHistsEkin(histseKin_top4);
				}
			} else (hevtsum->Fill(EVT_OTHER));
		}
	}
	if (pass) {
		EpProcessor::handle(d);
	}
}

void ProcTops::write(){
	if (hists!=NULL) {
		dirout->cd();
		hevtsum->Write();
		hists->Write();
		histseKin->Write();
		dirout->cd("top0");
		yields_top0->Write();
		histsMM_top0->Write();
		histseKin_top0->Write();
		dirout->cd("top1");
		yields_top1->Write();
		histsMM_top1->Write();
		histseKin_top1->Write();
		dirout->cd("top2");
		yields_top2->Write();
		histsMM_top2->Write();
		histseKin_top2->Write();
		dirout->cd("top3");
		yields_top3->Write();
		histsMM_top3->Write();
		histseKin_top3->Write();
		dirout->cd("top4");
		yields_top4->Write();
		histsMM_top4->Write();
		histseKin_top4->Write();
		//if (d->h10.is_sim) {
			//dirout->cd();
			if(dirout->cd("mc")){
				yields_mc->Write();
				histsMM_mc->Write();
				histseKin_mc->Write();
		    }
		//}
	}
}


void ProcTops::calcMM(DataOmega *d, Bool_t ismc /* = kFALSE */) {
	Data2pi *tp = &(d->twoPi);
	if (ismc) tp = &(d->twoPi_mc);
	tp->mm2ppippim = lvMMtop1.Mag2();
	tp->mm2ppip    = lvMMtop2.Mag2();
	tp->mm2ppim    = lvMMtop3.Mag2();
	tp->mm2pippim  = lvMMtop4.Mag2();
}

void ProcTops::calcVars(DataOmega* d, Bool_t ismc /* = kFALSE */){
	Data2pi *tp = &(d->twoPi);
	if (ismc) tp = &(d->twoPi_mc);
	
	lvQCMS   = lvQ;
	lvP0CMS  = _4vP0;
	lvPCMS   = lvP;
	lvPipCMS = lvPip;
	lvPimCMS = lvPim;
	lvQCMS.Boost(-1*lvW.BoostVector());
	lvP0CMS.Boost(-1*lvW.BoostVector());
	lvPCMS.Boost(-1*lvW.BoostVector());
	lvPipCMS.Boost(-1*lvW.BoostVector());
	lvPimCMS.Boost(-1*lvW.BoostVector());
	
	Float_t Mppip = (lvPCMS + lvPipCMS).Mag();
	Float_t Mppim = (lvPCMS + lvPimCMS).Mag();
	Float_t Mpippim = (lvPipCMS + lvPimCMS).Mag();
	
	tp->asgnmt1.M1 = Mppip;
	tp->asgnmt1.M2 = Mpippim;
	tp->asgnmt1.theta = getTheta(lvPimCMS);
	tp->asgnmt1.phi = getPhi(lvPimCMS);
	//tp->asgnmt1.alpha = getAlpha(lvPimCMS);
	tp->asgnmt1.alpha = 180;
	
	
	tp->asgnmt2.M1 = Mppip;
	tp->asgnmt2.M2 = Mpippim;
	tp->asgnmt2.theta = getTheta(lvPCMS);
	tp->asgnmt2.phi = getPhi(lvPCMS);
	//tp->asgnmt2.alpha = getAlpha(lvPCMS);
	tp->asgnmt2.alpha = 180;
	
	tp->asgnmt3.M1 = Mppip;
	tp->asgnmt3.M2 = Mppim;
	tp->asgnmt3.theta = getTheta(lvPipCMS);
	tp->asgnmt3.phi = getPhi(lvPipCMS);
	//tp->asgnmt1.alpha = getAlpha(lvPipCMS);
	tp->asgnmt3.alpha = 180;
	
	//helicity
	if(!ismc) {tp->h = d->h10.evthel;} //e1f
	
}

Float_t ProcTops::getTheta(TLorentzVector lv){
	Float_t retVal = 0;
	/* ***atrivedi 09-13-12*** */
	TVector3 lv_3vec = lv.Vect();
	TVector3 lvQCMS_3vec = lvQCMS.Vect();
	//retVal = RadToDeg()*ACos( (lv.Dot(lvQCMS))/(lv.P()*lvQCMS.P()) );
	retVal = RadToDeg()*ACos( (lv_3vec.Dot(lvQCMS_3vec))/(lv_3vec.Mag()*lvQCMS_3vec.Mag()) );
	/***************************/
	return retVal;
}

Float_t ProcTops::getPhi(TLorentzVector lv) {
	Float_t retVal = 0;
	retVal = RadToDeg()*invTan(lv.Py(), lv.Px());
	return retVal;
}

Float_t ProcTops::invTan(Float_t y, Float_t x){
	Float_t retVal = 0; 
	if (x > 0 && y > 0)  retVal = ATan(y/x);          //1st Quad.
	if (x < 0 && y > 0)  retVal = ATan(y/x) + Pi();   //2nd Quad
	if (x < 0 && y < 0)  retVal = ATan(y/x) + Pi();   //3rd Quad
	if (x > 0 && y < 0)  retVal = ATan(y/x) + 2*Pi(); //4th Quad 
	if (x == 0 && y > 0) retVal = Pi()/2;
	if (x == 0 && y < 0) retVal = 3*Pi()/2; 
	return retVal;  
}

void ProcTops::McKin(DataOmega *d) {
	lvQ.SetXYZT(0,0,0,0);
	lvW.SetXYZT(0,0,0,0);
	lvE.SetXYZT(0,0,0,0);
	lvP.SetXYZT(0,0,0,0);
	lvPip.SetXYZT(0,0,0,0);
	lvPim.SetXYZT(0,0,0,0);
	//printf("num mc = %i\n",d->h10.mcnentr);
	for (Int_t idx = 0; idx < d->h10.mcnentr; idx++) {
		Int_t _id = d->h10.mcid[idx];
		Float_t _p = d->h10.mcp[idx];
		Float_t _theta = d->h10.mctheta[idx]*DegToRad();
		Float_t _phi = d->h10.mcphi[idx]*DegToRad();
		Float_t _mass = d->h10.mcm[idx];
		Float_t _energy = Sqrt(_p*_p+_mass*_mass);
		Float_t _pz = _p*Cos(_theta);
		Float_t _py = _p*Sin(_phi)*Sin(_theta);
		Float_t _px = _p*Cos(_phi)*Sin(_theta);
		Int_t _sector = (d->h10.mcphi[idx]+30)/60.0+1;
		while (_sector<0) _sector+=6;
		while (_sector>6) _sector-=6;

		Float_t phiD = d->h10.mcphi[idx];
		Float_t thetaD = d->h10.mctheta[idx];
		TLorentzVector lvtmp;
		TLorentzVector lvtmp2;
		//printf("id = %i\n");
		switch(_id) {
		case ELECTRON:
			lvE.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PROTON:
			lvP.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PIP:
			lvPip.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PIM:
			lvPim.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		default:
			break;
		}
	}

	lvQ = _4vE0-lvE;
	lvW = lvQ+_4vP0;
	d->twoPi_mc.Q2 = -1*(lvQ.Mag2());
	d->twoPi_mc.W = lvW.Mag();
	
	d->updateEkin(kTRUE);
	calcMM(d,kTRUE);
	calcVars(d, kTRUE);
}

void ProcTops::set4vecHadrons(DataOmega* d){
	//Track Identified as hadrons
	if (d->h10idxP>0) {
		Double_t mom = d->h10.p[d->h10idxP];
		Double_t px = mom*d->h10.cx[d->h10idxP];
		Double_t py = mom*d->h10.cy[d->h10idxP];
		Double_t pz = mom*d->h10.cz[d->h10idxP];
		Double_t energy = Sqrt(mom*mom+MASS_P*MASS_P);
		lvP.SetPxPyPzE(px,py,pz,energy);
	}
	if (d->h10idxPip>0) {
		Double_t mom = d->h10.p[d->h10idxPip];
		Double_t px = mom*d->h10.cx[d->h10idxPip];
		Double_t py = mom*d->h10.cy[d->h10idxPip];
		Double_t pz = mom*d->h10.cz[d->h10idxPip];
		Double_t energy = Sqrt(mom*mom+MASS_PIP*MASS_PIP);
		lvPip.SetPxPyPzE(px,py,pz,energy);
	}
	if (d->h10idxPim>0) {
		Double_t mom = d->h10.p[d->h10idxPim];
		Double_t px = mom*d->h10.cx[d->h10idxPim];
		Double_t py = mom*d->h10.cy[d->h10idxPim];
		Double_t pz = mom*d->h10.cz[d->h10idxPim];
		Double_t energy = Sqrt(mom*mom+MASS_PIM*MASS_PIM);
		lvPim.SetPxPyPzE(px,py,pz,energy);
	}
	
	lvMMtop1 = (lvW-(lvP+lvPip+lvPim));
	lvMMtop2 = (lvW-(lvP+lvPip));
    lvMMtop3 = (lvW-(lvP+lvPim));  
	lvMMtop4 = (lvW-(lvPip+lvPim)); 
}
