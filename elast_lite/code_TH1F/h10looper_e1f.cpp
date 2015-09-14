#define h10looper_e1f_cxx
#include "h10looper_e1f.h"
#include "cuts.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

h10looper_e1f::h10looper_e1f(TString h10type,TChain* h10chain,TString fout_name,Long64_t nentries) : fChain(0) 
{
	Info("h10looper_e1f::h10looper_e1f","Setting up h10looper_e1f...\n");

	//! Decode h10type (proc_h10_lite already makes sure that the h10type is OK)
	TObjArray *h10type_tokens = h10type.Tokenize(":");
	_expt = h10type_tokens->At(0)->GetName();
	_dtyp = h10type_tokens->At(1)->GetName();
	_seq  = h10type_tokens->At(2)->GetName();

	//! Setup h10chain
	Init(h10chain);

	//! _fout
	_fout=new TFile(fout_name.Data(),"RECREATE");

	//! _nentries_to_proc
	Long64_t nentries_chain=h10chain->GetEntries();
	nentries==0?_nentries_to_proc=nentries_chain:_nentries_to_proc=nentries;
	Info("h10looper_e1f::h10looper_e1f", "Number of entries in specified by user =  %llu",nentries);
	Info("h10looper_e1f::h10looper_e1f", "Number of entries in TChain =  %llu",nentries_chain);
	Info("h10looper_e1f::h10looper_e1f", "Number of entries to processess =  %llu",_nentries_to_proc);

	//! set up eid
	Info("h10looper_e1f::h10looper_e1f", "Setting up eid pars for dtyp=%s",_dtyp.Data());
	setup_eid_cutpars(_dtyp);

	//! set up _pcorr
	Info("h10looper_e1f::h10looper_e1f", "Setting up _pcorr for dtyp=%s",_dtyp.Data());
	TString ws=getenv("WORKSPACE");;
	TString fpcorr=TString::Format("%s/ana2pi/MomCorr",ws.Data());
	_pcorr=new MomCorr_e1f((char *)fpcorr.Data());


	//!output objects
	//!event level
	_hevt=new TH1F("hevt","Event statistics",NUM_EVT_STATS,0.5,NUM_EVT_STATS+0.5);
	_hevt->GetXaxis()->SetBinLabel(EVT_TRG,"trgr");
	_hevt->GetXaxis()->SetBinLabel(EVT_E,"e^{-}");
	_hevt->GetXaxis()->SetBinLabel(EVT_E_INFID,"e^{-} infid");
	_hevt->GetXaxis()->SetBinLabel(EVT_ELSTC,"elstc-inc");
	if (_seq=="recon"){
		//! EID
		_heid=new TH1F("heid","EID statistics",NUM_EID_STATS,0.5,NUM_EID_STATS+0.5);
		_heid->GetXaxis()->SetBinLabel(EID_TRG,"trgr");
		_heid->GetXaxis()->SetBinLabel(EID_GPART0,"gpart>0");
		_heid->GetXaxis()->SetBinLabel(EID_Q,"q==-1");
		_heid->GetXaxis()->SetBinLabel(EID_HIT_DC,"hit_dc");
		_heid->GetXaxis()->SetBinLabel(EID_HIT_CC,"hit_cc");
		_heid->GetXaxis()->SetBinLabel(EID_HIT_SC,"hit_sc");
		_heid->GetXaxis()->SetBinLabel(EID_HIT_EC,"hit_ec");
		_heid->GetXaxis()->SetBinLabel(EID_STAT, "stat>0");
		_heid->GetXaxis()->SetBinLabel(EID_DC_STAT, "dc_stat>0");
		_heid->GetXaxis()->SetBinLabel(EID_P_MIN_ECTH,"p>p_ecth");
		_heid->GetXaxis()->SetBinLabel(EID_SF,"pass sf");
		//! EFID
		_hefid=new TH1F("hefid","EFID statistics",NUM_EFID_STATS,0.5,NUM_EFID_STATS+0.5);
		_hefid->GetXaxis()->SetBinLabel(EFID_TOT,"total");
		_hefid->GetXaxis()->SetBinLabel(EFID_IN,"infid");
		//! pcorr
		if (_dtyp=="exp"){
			//!pcorr
			_hpcorr_dpVp=new TH2F("hdpVp","#Deltap vs. p",550,0,5.5,160,-0.08,0.08);
			_hpcorr_dcx=new TH1F("hdcx", "#Deltacx",60,-0.01,0.01);
			_hpcorr_dcy=new TH1F("hdcy", "#Deltacy",60,-0.01,0.01);
			_hpcorr_dcz=new TH1F("hdcz", "#Deltacz",60,-0.01,0.01);
			_hpcorr_dp=new TH1F("hdp", "#Deltap", 160,-0.08,0.08);
		}
	}
	//! delast
	_hW=new TH1F("hW","hW",250,0,5);
	_helast=new TH1F("helast","helast",100,0.7,1.2);
	_hf=new TH1F*[6];
	_hc=new TH1F*[6];
	for (int isctr=0;isctr<6;isctr++){
		TString name;
		TString title;
		//! _hf
		name=TString::Format("hf_s%d",isctr+1);
		title=TString::Format("#phi=[%d,%d]",PHI_FULL_SECTOR[isctr][0],PHI_FULL_SECTOR[isctr][1]);
		_hf[isctr]=new TH1F(name,title,_NBINS,_THETA_MIN,_THETA_MAX);
		//! _hc
		name=TString::Format("hc_s%d",isctr+1);
		title=TString::Format("#phi=[%d,%d]",PHI_CENTRAL_SECTOR[isctr][0],PHI_CENTRAL_SECTOR[isctr][1]);
		_hc[isctr]=new TH1F(name,title,_NBINS,_THETA_MIN,_THETA_MAX);
	}
	/*printf("debug = %s\n",_hf[-1]->GetName());
	printf("debug = %s\n",_hf[-5]->GetName());*/
   
	//! Kinematics: initial and final: lvE0, lvP0
	_lvE0.SetXYZM(0,0,E0_P,MASS_E);
	_lvP0.SetXYZM(0,0,0,MASS_P);
	_lvE1=new TLorentzVector(); 
	_lvQ=new TLorentzVector();
	_lvW=new TLorentzVector(); 
	_Q2=_W=0;
	_theta=_phi=0;

	Info("h10looper_e1f::h10looper_e1f","Following information was setup:");
	Info("","---------------------------------------");
	Info("","expt:dtyp:seq=%s:%s:%s",_expt.Data(),_dtyp.Data(),_seq.Data());
	Info("","fout=%s",_fout->GetName());
	Info("","nentries_to_proc=%llu",_nentries_to_proc);
	Info("","Beam momentum and energy=%.3f,%.3f",_lvE0.Pz(),_lvE0.E());
	Info("","*** eid-cut pars ***");
	Info("","p_min_ECth=%.3f",_p_min_ECth);
	for (int isctr=0;isctr<6;isctr++){
		Info("","sf_low_sctr%d_pol3: %.5f,%.5f,%.5f,%.5f",
			isctr+1,_sf_min[isctr]->GetParameter(0),_sf_min[isctr]->GetParameter(1),_sf_min[isctr]->GetParameter(2),_sf_min[isctr]->GetParameter(3));
		Info("","sf_max_sctr%d_pol3: %.5f,%.5f,%.5f,%.5f",
			isctr+1,_sf_max[isctr]->GetParameter(0),_sf_max[isctr]->GetParameter(1),_sf_max[isctr]->GetParameter(2),_sf_max[isctr]->GetParameter(3));
	}
	Info("","---------------------------------------");
	Info("h10looper_e1f::h10looper_e1f","Done setting up h10looper_e1f\n");
	
}

h10looper_e1f::~h10looper_e1f()
{
	Info("h10looper_e1f::~h10looper_e1f","");
	if (!fChain) return;
	delete fChain->GetCurrentFile();

	delete[] _sf_mean;
	delete[] _sf_min;
	delete[] _sf_max;
	
	delete _hevt;
	delete _heid;
	delete _hefid;
	delete _hpcorr_dpVp;
	delete _hpcorr_dcx;
   delete _hpcorr_dcy;
   delete _hpcorr_dcz;
   delete _hpcorr_dp;
   delete _hW;
   delete _helast;
   delete[] _hf;
   delete[] _hc;

   delete _fout;

   delete _lvE1;
   delete _lvQ;
   delete _lvW;
}

void h10looper_e1f::setup_eid_cutpars(TString dtyp)
{
	_p_min_ECth=P_MIN_ECTH;
	_sf_mean=new TF1*[6];
	_sf_min=new TF1*[6];
	_sf_max=new TF1*[6];
	for (int isctr=0;isctr<6;isctr++){
		_sf_mean[isctr]=new TF1(TString::Format("sf_mean_s%d",isctr+1),"pol3",0,10);
		_sf_min[isctr]= new TF1(TString::Format("sf_min_s%d",isctr+1), "pol3",0,10);
		_sf_max[isctr]= new TF1(TString::Format("sf_max_s%d",isctr+1), "pol3",0,10);
		for (int ipar=0;ipar<4;ipar++){
			if (dtyp=="exp"){
				_sf_mean[isctr]->FixParameter(ipar,SF_MEAN_EXP[isctr][ipar]);
				_sf_min[isctr]->FixParameter(ipar,_sf_mean[isctr]->GetParameter(ipar)-3*SF_SIGMA_EXP[isctr][ipar]);
				_sf_max[isctr]->FixParameter(ipar,_sf_mean[isctr]->GetParameter(ipar)+3*SF_SIGMA_EXP[isctr][ipar]);
			}else if (dtyp=="sim"){
				_sf_mean[isctr]->FixParameter(ipar,SF_MEAN_SIM[isctr][ipar]);
				_sf_min[isctr]->FixParameter(ipar,_sf_mean[isctr]->GetParameter(ipar)-3*SF_SIGMA_SIM[isctr][ipar]);
				_sf_max[isctr]->FixParameter(ipar,_sf_mean[isctr]->GetParameter(ipar)+3*SF_SIGMA_SIM[isctr][ipar]);
			}
		}
	}
}

void h10looper_e1f::Loop(){
	Info("h10looper_e1f::Loop()","");

	if (fChain == 0) return;

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<_nentries_to_proc;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		if ( (jentry+1)%100000==0 ) Info("h10looper_e1f::Loop()","Processing event %llu",jentry+1);
		//! Begin processing event
		reset_ekin();
		if (_seq=="recon"){
			_hevt->Fill(EVT_TRG);
			if (evt_trigger_electron()){
				_hevt->Fill(EVT_E);
				set_ekin();
				if (electron_infid()){
					_hevt->Fill(EVT_E_INFID);
					if (_dtyp=="exp"){
						mom_corr_electron();
						set_ekin();
					}
					if (_W>_W_CUT_MIN && _W<_W_CUT_MAX){
						_hevt->Fill(EVT_ELSTC);
						make_delast();
					}
				}
			}
		}else if (_seq=="thrown"){
			_hevt->Fill(EVT_TRG);
			set_ekin();
			if (_W>_W_CUT_MIN && _W<_W_CUT_MAX){
				_hevt->Fill(EVT_ELSTC);
				make_delast();
			}
		}
	}//end event loop
	//! Print some statistics
	Info("h10looper_e1f::Loop()","*** statistics *** ");
	Info("","*** Total number of Elastic events from hevt = %.2f***",_hevt->GetBinContent(EVT_ELSTC));
	//! Now sum up Elastic events from individual sectors
	int nelast=0;
	for (int isctr=0;isctr<6;isctr++){
		nelast+=_hf[isctr]->GetEntries();
	}
	Info("","*** Total number of Elastic events after sum_X(hf_sX) = %.2f***",nelast);
	Info("h10looper_e1f::Loop()","****** ");
	
	//! Write _fout
	_fout->Write();
}//end Loop()

bool h10looper_e1f::evt_trigger_electron(){
	bool ret=kFALSE;

	int gprt=gpart;
	int chrg=q[0];
	bool hitDC=dc[0]>0,hitCC=cc[0]>0,hitSC=sc[0]>0,hitEC=ec[0]>0;
	int idxDC=dc[0]-1,idxCC=cc[0]-1,idxSC=sc[0]-1,idxEC=ec[0]-1;
	int stt=stat[0], dc_stt=dc_stat[idxDC];

	int sctr=ec_sect[idxEC];
	float mom=p[0];
	float sf=etot[idxEC]/mom;

	_heid->Fill(EID_TRG);
	bool pass_lvl_1=kFALSE;
	if (gprt>0){
		_heid->Fill(EID_GPART0);
		if (chrg==-1){
			_heid->Fill(EID_Q);
			if (hitDC){
				_heid->Fill(EID_HIT_DC);
				if (hitCC){
					_heid->Fill(EID_HIT_CC);
					if (hitSC){
						_heid->Fill(EID_HIT_SC);
						if (hitEC){
							_heid->Fill(EID_HIT_EC);
							if (stt>0){
								_heid->Fill(EID_STAT);
								if (dc_stt>0){
									_heid->Fill(EID_DC_STAT);
									pass_lvl_1=kTRUE;
								}
							}
						}
					}
				}
			}
		}
	} 

	if (pass_lvl_1){
		if (mom>_p_min_ECth){
			_heid->Fill(EID_P_MIN_ECTH);
			float sf_min=_sf_min[sctr-1]->Eval(mom);
			float sf_max=_sf_max[sctr-1]->Eval(mom);
			if (sf>sf_min && sf<sf_max){
				_heid->Fill(EID_SF);
				ret=kTRUE;
			}
		}
	}

	return ret;
}

bool h10looper_e1f::electron_infid(){
	Int_t id=ELECTRON;
	Double_t mom=p[0];
	Int_t ecidx=ec[0]-1;
	Int_t sector=ec_sect[ecidx];
	Int_t paddle = 0;
	
	_hefid->Fill(EFID_TOT);
	bool ret=Cuts::Fiducial(id, mom, _theta, _phi, sector, paddle);
	if (ret==kTRUE){
		_hefid->Fill(EFID_IN);
	}
	return ret;
}

void h10looper_e1f::mom_corr_electron(){
	//! First store uncorr mom
	float p_uncorr=_lvE1->P();
	float cx_uncorr=_lvE1->Px()/_lvE1->P();
	float cy_uncorr=_lvE1->Py()/_lvE1->P();
	float cz_uncorr=_lvE1->Pz()/_lvE1->P();
	//! Now correct
	int q=-1;
	int id=ELECTRON;
	TLorentzVector pcorr=_pcorr->PcorN(*_lvE1,q,id);
	//! update h10 data
	p[0]=pcorr.P();
	cx[0]=pcorr.Px()/pcorr.P();
	cy[0]=pcorr.Py()/pcorr.P();
	cz[0]=pcorr.Pz()/pcorr.P();
	//! update pcorr objects
	_hpcorr_dpVp->Fill(p[0],p[0]-p_uncorr);
	_hpcorr_dcx->Fill(cx[0]-cx_uncorr);
	_hpcorr_dcy->Fill(cy[0]-cy_uncorr);
	_hpcorr_dcz->Fill(cz[0]-cz_uncorr);
	_hpcorr_dp->Fill(p[0]-p_uncorr);
}

void h10looper_e1f::set_ekin(){
	if (_seq=="recon"){
		float px=p[0]*cx[0];
		float py=p[0]*cy[0];
		float pz=p[0]*cz[0];
		_lvE1->SetXYZM(px,py,pz,MASS_E);
	}else if (_seq=="thrown"){
		for (int inprt=0; inprt<nprt;inprt++){
			if (pidpart[inprt]==3){
				float px=pxpart[inprt];
				float py=pypart[inprt];
				float pz=pzpart[inprt];
				float e=epart[inprt];
				_lvE1->SetPxPyPzE(px,py,pz,e);
			}
		}
	}
	TLorentzVector lvQ=_lvE0-(*_lvE1);
	TLorentzVector lvW=lvQ+_lvP0;
	*_lvQ=lvQ;
	*_lvW=lvW;
	_Q2=-1*_lvQ->Mag2();
	_W=_lvW->Mag();
	_theta=_lvE1->Theta()*TMath::RadToDeg();
	float phi=_lvE1->Phi()*TMath::RadToDeg();// [-180,180]
	_phi=phi<-30?phi+360:phi; // [-30,330]
}

void h10looper_e1f::reset_ekin(){
	_lvE1->SetXYZT(0,0,0,0);
	_lvQ->SetXYZT(0,0,0,0);
	_lvW->SetXYZT(0,0,0,0);
	_Q2=_W=0;
	_theta=_phi=0;
}

void h10looper_e1f::make_delast(){
	//! _hW and _helastic
	_hW->Fill(_W);
	_helast->Fill(_W);
	int isctr=get_sector()-1;
	//Info("h10looper_e1f::make_delast","isctr=%d",isctr);
	//! _hf
	_hf[isctr]->Fill(_theta);
	//! _hc
	if (_phi>=PHI_CENTRAL_SECTOR[isctr][0] && _phi<=PHI_CENTRAL_SECTOR[isctr][1]){
		_hc[isctr]->Fill(_theta);
	}
	return;
}

int h10looper_e1f::get_sector(){
	int sector=-1;
	if ( (_phi>=PHI_FULL_SECTOR[0][0]  && _phi<PHI_FULL_SECTOR[0][1]) || _phi==PHI_FULL_SECTOR[5][1]){
		sector=1;
	}else if (_phi>=PHI_FULL_SECTOR[1][0] && _phi<PHI_FULL_SECTOR[1][1]){
		sector=2;
	}else if (_phi>=PHI_FULL_SECTOR[2][0] && _phi<PHI_FULL_SECTOR[2][1]){
		sector=3;
	}else if (_phi>=PHI_FULL_SECTOR[3][0] && _phi<PHI_FULL_SECTOR[3][1]){
		sector=4;
	}else if (_phi>=PHI_FULL_SECTOR[4][0] && _phi<PHI_FULL_SECTOR[4][1]){
		sector=5;
	}else if (_phi>=PHI_FULL_SECTOR[5][0] && _phi<PHI_FULL_SECTOR[5][1]){
		sector=6;
	}else{
		Info("h10looper_e1f::get_sector","sector not found for electron phi,theta,Q2,W=%.2f,%.2f,%.2f,%.2f. Returning sector=-1",_phi,_theta,_Q2,_W);
		Info("h10looper_e1f::get_sector","file=%s\n",fChain->GetCurrentFile()->GetName());
		/*printf("nprt=%d\n",nprt);
		for (int inprt=0; inprt<nprt;inprt++){
			printf("pidpart[%d]=%d\n",inprt,pidpart[inprt]);
			printf("px[%d]=%.2f\n",inprt,pxpart[inprt]);
			printf("py[%d]=%.2f\n",inprt,pypart[inprt]);
			printf("pz[%d]=%.2f\n",inprt,pzpart[inprt]);
			printf("e[%d]=%.2f\n",inprt,epart[inprt]);
		}*/
	}
	return sector;
}
