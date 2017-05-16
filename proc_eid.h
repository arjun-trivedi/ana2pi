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
#include <stdio.h>
#include "wrpr_vertex_e16.h" //! [01-15-16] for zvtx corr

/*
[01-17-16]
+ Added corr_zvtx() which is currently on setup to work for e16:ER and therefore
  is called in handle() only if h10->expt=="e16" and h10->dtyp=="exp" 
*/

 /*
[07-08-16]
+ Added relevant functions and objects, FOR E16 ONLY, to make zvtx and ECin_min cut.
+ Added objects to monitor zvtx corr (which was added on 01-17-16)
*/

/*
[07-15-16]
+ Added options to not apply cuts that need to be studied.
+ Currently the only 3 cuts (*nphe is excluded for now*) that need to be, and have been, studied are ECfid, zvtx, and SF.
	+ 1. ECfid: added user option '_study_ECfid' (Note if applied ECfid-cut=E1F values)
	+ 2. zvtx:  added used option '_study_zvtx'
	+ 3. SF:    this already has its own specialized processor 'study_eid_play.h', which was made earlier, before I though
	            of this new way.

+ [01-15-17] _eidTool is now updated to support all e16 cuts. In fact, now it is e1f that is not
  	         optimal, but that is OK because I am concentration only on e16 for now.
*/

using namespace TMath;
using namespace ParticleConstants;
using namespace AnalysisConstants;

class ProcEid : public EpProcessor
{
public:
	ProcEid(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		    Bool_t make_tree=kFALSE, bool study_ECfid=kFALSE,bool study_zvtx=kFALSE, bool study_nphe=kFALSE);
	ProcEid(DataH10* dataH10,DataAna* dataAna);
	~ProcEid();
	
	void handle();
	//void write();
	Bool_t goodE();
	Bool_t goodE_bos();
	void updateEid();
	void updateEid_zvtx_only();
	void updateEkin(Bool_t useMc=kFALSE,Bool_t McHasPARTBanks=kFALSE);
	Float_t getCCtheta(Float_t x_sc, Float_t y_sc, Float_t z_sc, Float_t cx_sc, Float_t cy_sc, Float_t cz_sc);
	void GetUVW(float xyz[3], float uvw[3]);
	void corr_zvtx();
	
protected:
	Eid* _eidTool;
	TDirectory* _dirmon;
	TDirectory* _dirzvtxcorr;
	TDirectory* _dircut;
	TTree* _t[NPROCMODES];
	//! _hzvtxcorr[6][2] for zvtxcorr monitoring, before and after cut
    TH1F*** _hzvtxcorr;
				
	static const Int_t NUM_EVTCUTS = 16;
	enum { EVT_NULL, EVT_TRIG, EVT_GPART1, EVT_STAT1, EVT_Q1,
	       EVT_DC1, EVT_CC1, EVT_SC1, EVT_EC1, EVT_NPHE,
	       EVT_DCSTAT1, EVT_ECLOW1, EVT_ECFID, EVT_ECIN_MIN, EVT_ZVTX, EVT_SF, EVT_BOS11
	     };
	bool _make_tree;
	bool _study_ECfid;
	bool _study_zvtx;
	bool _study_nphe;
};

ProcEid::ProcEid(TDirectory *td, DataH10* dataH10, DataAna* dataAna, 
                 Bool_t make_tree/* = kFALSE*/,bool study_ECfid/*=kFALSE*/,bool study_zvtx/*=kFALSE*/,bool study_nphe/*=kFALSE*/)
                 :EpProcessor(td, dataH10, dataAna)
{
	TString path;
  	path=getenv("WORKSPACE");
  	
    //![11-17-15] Use same eid-cut pars for E16 and E1F
    //![02-17-16]
    //! + same eid-cut pars for sim:{e16,e1f}
   	//! + For exp:e16, I created eid.exp.e16.out where I have updated the SFpars from 
    //!   'study_eid/study_SF/results_SFvp_e16/cutpars/exp_peakSF.txt' 
    //!   (Used study_eid/study_SF/dev/SFpars_H-L_2_MU-SG.py to convert pars in above file from cut_h/l to mean/sgma)
    /*if      (dH10->dtyp=="sim") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.mc.out",path.Data())).Data());
	else if (dH10->dtyp=="exp") {
		if      (dH10->expt=="e1f") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.exp.out",path.Data())).Data());
		else if (dH10->expt=="e16") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.exp.e16.out",path.Data())).Data());
	}*/
	//! 01-15-17] Initialize _eidTool in *dependence* of expt (For details see comments on top of file)
  	if (dH10->expt=="e1f") {
		if      (dH10->dtyp=="sim") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.mc.out",path.Data())).Data());
		else if (dH10->dtyp=="exp") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.exp.out",path.Data())).Data());
	}
	else if (dH10->expt=="e16") {
		if      (dH10->dtyp=="sim") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.mc.e16.out",path.Data())).Data());
		else if (dH10->dtyp=="exp") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.exp.e16.out",path.Data())).Data());
	}else  Info("ProcEid::ProcEid()", "_eidTool not initialized");
	
	_make_tree=make_tree;
	_study_ECfid=study_ECfid;
	_study_zvtx=study_zvtx;
	_study_nphe=study_nphe;

	if (_study_ECfid){
		Info("ProcEid::ProcEid()", "_study_ECfid=kTRUE, therefore this cut will not be applied");
	}
	if (_study_zvtx){
		Info("ProcEid::ProcEid()", "_study_zvtx=kTRUE, therefore this cut will not be applied");
	}
	if (_study_nphe){
		Info("ProcEid::ProcEid()", "_study_nphe=kTRUE, therefore this cut will not be applied");
	}


	_dirmon=NULL;
	_dirzvtxcorr=NULL;
	_dircut=NULL;

	td->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	//hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT_TRIG,"Trigger");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPART1,"gpart>0");
	hevtsum->GetXaxis()->SetBinLabel(EVT_STAT1,"stat>0");
	hevtsum->GetXaxis()->SetBinLabel(EVT_Q1,"q=-1");
	hevtsum->GetXaxis()->SetBinLabel(EVT_DC1,"DC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_CC1,"CC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_SC1,"SC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_EC1,"EC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_NPHE,"NPHE");
	hevtsum->GetXaxis()->SetBinLabel(EVT_DCSTAT1,"dc_stat>0");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ECLOW1,"EC Threshold");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ECFID,"EC Fid.");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ECIN_MIN,"ec_ei > ECmin");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ZVTX,"pass zvrtx");
	hevtsum->GetXaxis()->SetBinLabel(EVT_SF,"SF");
	hevtsum->GetXaxis()->SetBinLabel(EVT_BOS11,"EVNT.id=11");
	hevtsum->SetMinimum(0.);

	//! Make monitor output objects
	_dirmon = dirout->mkdir(TString::Format("monitor"));
	dAna->makeHistsEid(hists[MONMODE][EVTINC], _dirmon);
	dAna->makeHistsEkin(histsEkin[MONMODE][EVTINC], _dirmon);

	//! 07-06-16 Make hists to monitor zvtxcorr
	_dirzvtxcorr=dirout->mkdir(TString::Format("zvtxcorr"));
	_dirzvtxcorr->cd();
	_hzvtxcorr=new TH1F**[6];
	for (int i=0;i<6;i++){
		_hzvtxcorr[i]=new TH1F*[2];
		for (int j=0;j<2;j++){
			TString name_sfx;
			if      (j==0) name_sfx="prec";
			else if (j==1) name_sfx="pstc";
			TString name=TString::Format("zvtxcorr_s%d_%s",i+1,name_sfx.Data());
			_hzvtxcorr[i][j]=new TH1F(name,name,100,-12,6);
		}
	}
	
	//! monmode tree	
	if (_make_tree){
		_dirmon->cd();
		_t[MONMODE] = new TTree("t","TTree containing data from ProdEid");
		dAna->addBranches_DataEkin(_t[MONMODE]);
		dAna->addBranches_DataEid(_t[MONMODE]);
	}

	//! Make cut output objects
	_dircut = dirout->mkdir(TString::Format("cut"));
	dAna->makeHistsEid(hists[CUTMODE][EVTINC], _dircut);
	dAna->makeHistsEkin(histsEkin[CUTMODE][EVTINC], _dircut);
	
	if (_make_tree){
		_dircut->cd();
		_t[CUTMODE] = new TTree("t","TTree containing data from ProdEid");
		dAna->addBranches_DataEkin(_t[CUTMODE]);
		dAna->addBranches_DataEid(_t[CUTMODE]);
	}

	//! [01-17-16] if e16:ER then corr_zvtx() 
    if (dH10->expt=="e16" and dH10->dtyp=="exp") {
    	Info("ProcEid::ProcEid()", "expt:exp=e16:ER and therefore will do corr_zvtx()");
    }
}

ProcEid::ProcEid(DataH10* dataH10, DataAna* dataAna)
                 :EpProcessor(dataH10, dataAna)
{
	TString path;
  	path=getenv("WORKSPACE");
  	
  	//![11-17-15] Use same eid-cut pars for E16 and E1F
    //![02-17-16]
    //! + same eid-cut pars for sim:{e16,e1f}
   	//! + For exp:e16, I created eid.exp.e16.out where I have updated the SFpars from 
    //!   'study_eid/study_SF/results_SFvp_e16/cutpars/exp_peakSF.txt' 
    //!   (Used study_eid/study_SF/dev/SFpars_H-L_2_MU-SG.py to convert pars in above file from cut_h/l to mean/sgma)
    /*if      (dH10->dtyp=="sim") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.mc.out",path.Data())).Data());
	else if (dH10->dtyp=="exp") {
		if      (dH10->expt=="e1f") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.exp.out",path.Data())).Data());
		else if (dH10->expt=="e16") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.exp.e16.out",path.Data())).Data());
	}*/
	//! 01-15-17] Initialize _eidTool in *dependence* of expt (For details see comments on top of file)
  	if (dH10->expt=="e1f") {
		if      (dH10->dtyp=="sim") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.mc.out",path.Data())).Data());
		else if (dH10->dtyp=="exp") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.exp.out",path.Data())).Data());
	}
	else if (dH10->expt=="e16") {
		if      (dH10->dtyp=="sim") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.mc.e16.out",path.Data())).Data());
		else if (dH10->dtyp=="exp") _eidTool = new Eid((char *)(TString::Format("%s/ana2pi/eid/eid.exp.e16.out",path.Data())).Data());
	}else  Info("ProcEid::ProcEid()", "_eidTool not initialized");

	//! [05-15-17] 
	//! + 'hevtsum' was earlier not created for this constructor
	//! + It is now added, but since 'td' does not exists for this case, 'eid' is made explicit in the 
	//! name of the histogram
	//td->cd();
	hevtsum = new TH1D("eid_hevtsum","Event Statistics from ProcEid",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	//hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT_TRIG,"Trigger");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPART1,"gpart>0");
	hevtsum->GetXaxis()->SetBinLabel(EVT_STAT1,"stat>0");
	hevtsum->GetXaxis()->SetBinLabel(EVT_Q1,"q=-1");
	hevtsum->GetXaxis()->SetBinLabel(EVT_DC1,"DC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_CC1,"CC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_SC1,"SC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_EC1,"EC");
	hevtsum->GetXaxis()->SetBinLabel(EVT_NPHE,"NPHE");
	hevtsum->GetXaxis()->SetBinLabel(EVT_DCSTAT1,"dc_stat>0");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ECLOW1,"EC Threshold");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ECFID,"EC Fid.");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ECIN_MIN,"ec_ei > ECmin");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ZVTX,"pass zvrtx");
	hevtsum->GetXaxis()->SetBinLabel(EVT_SF,"SF");
	hevtsum->GetXaxis()->SetBinLabel(EVT_BOS11,"EVNT.id=11");
	hevtsum->SetMinimum(0.);
}

ProcEid::~ProcEid(){
	delete _eidTool;
	delete _t;
	delete _dirmon;
	delete _dirzvtxcorr;
	delete _dircut;
	delete[] _hzvtxcorr;
}
	
void ProcEid::handle() {
	//Info("ProcEid::handle()", "");
	pass = kFALSE;
	
	hevtsum->Fill(EVT_TRIG);
	
	//! Monitor
	updateEid();
	updateEkin();
	dAna->fillHistsEid(hists[MONMODE][EVTINC]);
	dAna->fillHistsEkin(histsEkin[MONMODE][EVTINC]);
	if (_make_tree){
		_t[MONMODE]->Fill();
	}
		
	//! Cut	
	Bool_t gE = kFALSE;

    /*if      (dH10->expt=="e1f" && _eidTool->eidParFileFound)  gE =  goodE();
    else if (dH10->expt=="e1f" && !_eidTool->eidParFileFound) gE =  goodE_bos();
    else if (dH10->expt=="e16" && _eidTool->eidParFileFound) gE =  goodE();
    else if (dH10->expt=="e16")                               gE =  goodE_bos(); //pars for e1-6 not yet obtained*/

    //! [01-17-16] if e16:ER then corr_zvtx() and update DataEid zvtx 
    if (dH10->expt=="e16" and dH10->dtyp=="exp") {
    	//! 1. Fill prec monitoring
    	int sctr=dAna->eid.sector;
		if (sctr!=0) _hzvtxcorr[sctr-1][0]->Fill(dAna->eid.vz);
		//! 2. do corr_zvtx()
    	corr_zvtx();
    	//! 3. update DataEid::vz
    	updateEid_zvtx_only();
    	//! 4. pstc monitoring
    	if (sctr!=0)_hzvtxcorr[sctr-1][1]->Fill(dAna->eid.vz);
	}

    gE=goodE();
    
	if (gE) {	
		dAna->fillHistsEid(hists[CUTMODE][EVTINC]);
		dAna->fillHistsEkin(histsEkin[CUTMODE][EVTINC]);

		if (_make_tree){
			_t[CUTMODE]->Fill();
		}
		
		dH10->id[0] = ELECTRON;
		pass = kTRUE;
		EpProcessor::handle();
	}
}

Bool_t ProcEid::goodE(){
	Bool_t retval = kFALSE;

	DataEid* eid = &dAna->eid;
	
	if (dH10->dtyp!="sim") { //atrivedi 020313 till _eidTool is fixed to use eid.mc.out
		if (dH10->id[0]==ELECTRON) hevtsum->Fill(EVT_BOS11);
	}
	if (dH10->gpart>0) {
		hevtsum->Fill(EVT_GPART1);
		if (dH10->stat[0]>0) {//(dH10->sc_sect[dH10->sc[0]-1]==dH10->ec_sect[dH10->ec[0]-1]){//if (dH10->stat[0]>0) {
			hevtsum->Fill(EVT_STAT1);
			if (dH10->q[0]==-1) {
				hevtsum->Fill(EVT_Q1);
				if (dH10->dc[0]>0) {
					hevtsum->Fill(EVT_DC1);
					if (dH10->cc[0]>0) { //(dH10->cc[0]>0 || dH10->dtyp=="sim") {
						hevtsum->Fill(EVT_CC1);
						if (dH10->sc[0]>0) {
							hevtsum->Fill(EVT_SC1);
							if (dH10->ec[0]>0) {
								hevtsum->Fill(EVT_EC1);
								if( (_study_nphe) || (!_study_nphe && (dH10->nphe[dH10->cc[0]-1]>20 || dH10->dtyp=="sim")) ){
									hevtsum->Fill(EVT_NPHE);
									if (dH10->dc_stat[dH10->dc[0]-1]>0) {
										hevtsum->Fill(EVT_DCSTAT1);
										if (_eidTool->PassThreshold(eid->p)) {
											hevtsum->Fill(EVT_ECLOW1);
											float uvw[3]={eid->ecU,eid->ecV,eid->ecW};
											if ( (_study_ECfid) || (!_study_ECfid && _eidTool->PassECFid(uvw)) ){
												hevtsum->Fill(EVT_ECFID);
												Int_t sector = eid->sector;
												if ( (dH10->expt!="e16") || (dH10->expt=="e16" && _eidTool->pass_ECin_min(sector,eid->ec_ei)) ){
													hevtsum->Fill(EVT_ECIN_MIN);
													if ( (dH10->expt!="e16" || _study_zvtx) || (dH10->expt=="e16" && !_study_zvtx && _eidTool->pass_zvtx(sector, eid->vz)) ){
														hevtsum->Fill(EVT_ZVTX);
														Float_t p = eid->p;
														Float_t sf = eid->etot/p;
														if (_eidTool->PassSF(sector,p,sf)) {
															hevtsum->Fill(EVT_SF);
															if (!_eidTool->Pass(sector,p,sf)) {
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
				}
			}
		}
	}
	return retval;
}

Bool_t ProcEid::goodE_bos(){
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

void ProcEid::updateEid(){
		
	//! From cooking
	dAna->eid.id = dH10->id[0];
	//! From DC
	dAna->eid.p = dH10->p[0];
	dAna->eid.dc_xsc = dH10->dc_xsc[dH10->dc[0]-1];
	dAna->eid.dc_ysc = dH10->dc_ysc[dH10->dc[0]-1];
	dAna->eid.dc_zsc = dH10->dc_zsc[dH10->dc[0]-1];
	dAna->eid.vz     = dH10->vz[0];
	//! From SC
	dAna->eid.sector = dH10->sc_sect[dH10->sc[0]-1];

	Float_t p = dH10->p[0];
	Float_t l_e = dH10->sc_r[dH10->sc[0]-1];
	Float_t t_e = dH10->sc_t[dH10->sc[0]-1];
	Float_t tOFF = t_e-(l_e/SOL);
	
	Float_t b = ( l_e/(t_e - tOFF) )/SOL; // = 1
	Float_t b_e = TMath::Sqrt((p*p)/(MASS_E*MASS_E+p*p));
	Float_t dt_e = l_e/(b_e*SOL) + tOFF - t_e; // = 0
	
	dAna->eid.b = b;
	dAna->eid.b_e = b_e;
	dAna->eid.dt_e = dt_e;
		
	//from EC
	dAna->eid.ec_ei = dH10->ec_ei[dH10->ec[0]-1];
	dAna->eid.ec_eo = dH10->ec_eo[dH10->ec[0]-1];
	dAna->eid.etot  = dH10->etot[dH10->ec[0]-1];
	Float_t ech_x = dH10->ech_x[dH10->ec[0]-1];
	Float_t ech_y = dH10->ech_y[dH10->ec[0]-1];
	Float_t ech_z = dH10->ech_z[dH10->ec[0]-1];
	dAna->eid.ech_x=ech_x;
	dAna->eid.ech_y=ech_y;
	dAna->eid.ech_z=ech_z;
	Float_t xyz[3]={ech_x,ech_y,ech_z};
	Float_t uvw[3]={0,0,0};
	GetUVW(xyz,uvw);
	dAna->eid.ecU=uvw[0];
	dAna->eid.ecV=uvw[1];
	dAna->eid.ecW=uvw[2];
	//from CC
	dAna->eid.cc      = dH10->cc[0];
	dAna->eid.nphe    = dH10->nphe[dH10->cc[0]-1];
	dAna->eid.cc_segm = (dH10->cc_segm[dH10->cc[0]-1]%1000)/10;
	dAna->eid.pmt     = (dH10->cc_segm[dH10->cc[0]-1]/1000)-1;
	//calculate cc_theta
	Float_t dc_xsc  = dH10->dc_xsc[dH10->dc[0]-1];
	Float_t dc_ysc  = dH10->dc_ysc[dH10->dc[0]-1];
	Float_t dc_zsc  = dH10->dc_zsc[dH10->dc[0]-1];
	Float_t dc_cxsc = dH10->dc_cxsc[dH10->dc[0]-1];
	Float_t dc_cysc = dH10->dc_cysc[dH10->dc[0]-1];
	Float_t dc_czsc = dH10->dc_czsc[dH10->dc[0]-1];
	dAna->eid.cc_theta = getCCtheta(dc_xsc, dc_ysc, dc_zsc, dc_cxsc, dc_cysc, dc_czsc);
}

void ProcEid::updateEid_zvtx_only(){
	dAna->eid.vz = dH10->vz[0];
}

void ProcEid::updateEkin(Bool_t useMc /*= kFALSE*/,Bool_t McHasPARTBanks/*=kFALSE*/) {
	const TLorentzVector _4vE0 = dH10->lvE0;
	const TLorentzVector _4vP0 = dH10->lvP0;
	TLorentzVector _4vE1;
		
	DataEkin *ekin = &dAna->eKin;
	if (useMc) ekin = &dAna->eKin_mc;
	
	if (!useMc){
		//! Kinematics
		Float_t mom = dH10->p[0];
		Float_t px = mom*dH10->cx[0];
		Float_t py = mom*dH10->cy[0];
		Float_t pz = mom*dH10->cz[0];
		_4vE1.SetPxPyPzE(px,py,pz,Sqrt(mom*mom+MASS_E*MASS_E));
		//! Vertex information
		ekin->vx=dH10->vx[0];
		ekin->vy=dH10->vy[0];
		ekin->vz=dH10->vz[0];
	} else{//! if useMc
		if (McHasPARTBanks){
			for (Int_t inprt=0; inprt<dH10->nprt;inprt++)	{
				if (dH10->pidpart[inprt]==3){
					Float_t px=dH10->pxpart[inprt];
					Float_t py=dH10->pypart[inprt];
					Float_t pz=dH10->pzpart[inprt];
					Float_t e=dH10->epart[inprt];
					_4vE1.SetPxPyPzE(px,py,pz,e);
					//! Vertex information
					ekin->vx=dH10->xpart[inprt];
					ekin->vy=dH10->ypart[inprt];
					ekin->vz=dH10->zpart[inprt];
				}
			}
		}else{//!McHasMCTKBans
			for (Int_t idx = 0; idx < dH10->mcnentr; idx++) {
				Int_t id = dH10->mcid[idx];
				if (id != ELECTRON) continue;
				Float_t mom = dH10->mcp[idx];
				Float_t theta = dH10->mctheta[idx]*DegToRad();
				Float_t phi = dH10->mcphi[idx]*DegToRad();
				Float_t px = mom*Cos(phi)*Sin(theta);
				Float_t py = mom*Sin(phi)*Sin(theta);
				Float_t pz = mom*Cos(theta);
				_4vE1.SetPxPyPzE(px,py,pz,Sqrt(mom*mom+MASS_E*MASS_E));
			}
			//! Vertex information (Note that for MCTK, only e- vertex information in MCTK banks)
			ekin->vx=dH10->mcvx_x_el;
			ekin->vy=dH10->mcvx_y_el;
			ekin->vz=dH10->mcvx_z_el;
		}
	}
	
	if (!useMc) ekin->sector=dH10->sc_sect[dH10->sc[0]-1];
	TLorentzVector _4vQ = _4vE0-_4vE1;
	if (!useMc) {ekin->sector = dH10->sc_sect[dH10->sc[0]-1];}
	ekin->W = (_4vQ+_4vP0).Mag();
	ekin->Q2 = -1*_4vQ.Mag2();
	ekin->nu = _4vQ.E();
	ekin->xb = ekin->Q2/(2*MASS_P*ekin->nu);
	ekin->E1 = _4vE1.E();
	ekin->theta1 = _4vE1.Theta()*RadToDeg();
	Double_t phitmp = _4vE1.Phi()*RadToDeg(); //phitmp = [-180,180]
	ekin->phi1 = phitmp<-30 ? phitmp+360 : phitmp; // [-30,330]
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

/*****************************************************
[07-09-15]
The following function, taken directly as sent to my by Evan,
takes the global x,y,z coordinates of the EC (ech_x,ech_y,ech_z)
and fills calculates the local U,V,W coordinates of the EC.
*****************************************************/
void ProcEid::GetUVW(float xyz[3], float uvw[3]) {
  enum { X, Y, Z };
  enum { U, V, W };
  Float_t phi, ec_phi, ec_the, tgrho, sinrho, cosrho;
  Float_t ylow, yhi, xi, yi, zi;
  ec_the = 0.4363323;
  ylow = -182.974;
  yhi = 189.956;
  tgrho = 1.95325;
  sinrho = 0.8901256;
  cosrho = 0.455715;
  //----------------
  phi = ATan2(xyz[Y], xyz[X]) * 57.29578;
  if (phi < 0.)
    phi += 360.;
  phi = phi + 30.;
  if (phi > 360.)
    phi -= 360.;
  ec_phi = 1.0471975 * int(phi / 60.);
  //----------------
  float rot11 = cos(ec_the) * cos(ec_phi);
  float rot12 = -sin(ec_phi);
  float rot13 = sin(ec_the) * cos(ec_phi);
  float rot21 = cos(ec_the) * sin(ec_phi);
  float rot22 = cos(ec_phi);
  float rot23 = sin(ec_the) * sin(ec_phi);
  float rot31 = -sin(ec_the);
  float rot32 = 0.;
  float rot33 = cos(ec_the);
  //----------------
  yi = xyz[X] * rot11 + xyz[Y] * rot21 + xyz[Z] * rot31;
  xi = xyz[X] * rot12 + xyz[Y] * rot22 + xyz[Z] * rot32;
  zi = xyz[X] * rot13 + xyz[Y] * rot23 + xyz[Z] * rot33;
  zi -= 510.32;
  //----------------
  uvw[U] = (yi - ylow) / sinrho;
  uvw[V] = (yhi - ylow) / tgrho - xi + (yhi - yi) / tgrho;
  uvw[W] = ((yhi-ylow)/tgrho+xi+(yhi-yi)/tgrho)/2./cosrho;
}

/*
[01-17-16] 
Currently only for e16:ER
*/
void ProcEid::corr_zvtx(){
	//! Get all required quanties 
	float p = dH10->p[0];
	float cx = dH10->cx[0];
	float cy = dH10->cy[0];
	float cz = dH10->cz[0];
	float px = p*cx;
	float py = p*cy;
	float pz = p*cz;
	float vx= dH10->vx[0];
	float vy= dH10->vy[0];
	float vz= dH10->vz[0];
	//! + Prepare to call .f routines
	float vx_corr=0, vy_corr=0,vz_corr=0;
	vertex_e16(px,py,pz,vx,vy,vz,vx_corr,vy_corr,vz_corr);
	//! Update dH10
	dH10->vx[0]=vx_corr;
	dH10->vy[0]=vy_corr;
	dH10->vz[0]=vz_corr;

	return;
}
#endif // PROCEID_H
