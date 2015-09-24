#define h10looper_e1f_cxx
#include "h10looper_e1f.h"
#include "cuts.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

h10looper_e1f::h10looper_e1f(TString h10type,TChain* h10chain,TString fout_name,Long64_t nentries,
						     TString adtnl_cut_opt) : fChain(0) 
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

	//! cuts to be made in addition to 'dflt'
    setup_adtnl_cut_opts(adtnl_cut_opt);
	
	//! set up eid
	Info("h10looper_e1f::h10looper_e1f", "Setting up eid pars for dtyp=%s",_dtyp.Data());
	setup_eid_cutpars(_dtyp);

	//! set up pid
	Info("h10looper_e1f::h10looper_e1f", "Setting up _pid_tool for dtyp=%s",_dtyp.Data());
	_pid_tool=new Pid(_dtyp);

	//! set up _pcorr
	Info("h10looper_e1f::h10looper_e1f", "Setting up _pcorr");
	TString ws=getenv("WORKSPACE");;
	TString fpcorr=TString::Format("%s/ana2pi/MomCorr",ws.Data());
	_pcorr=new MomCorr_e1f((char *)fpcorr.Data());


	//!output objects
	//!event level
	_hevt=new TH1D("hevt","Event statistics",NUM_EVT_STATS,0.5,NUM_EVT_STATS+0.5);
	_hevt->GetXaxis()->SetBinLabel(EVT_TRG,"trgr");
	_hevt->GetXaxis()->SetBinLabel(EVT_E,"e^{-}");
	_hevt->GetXaxis()->SetBinLabel(EVT_E_INFID,"e^{-} infid");
	_hevt->GetXaxis()->SetBinLabel(EVT_ELSTC,"elstc-inc");
	if (_seq=="recon"){
		//! EID
		_heid=new TH1D("heid","EID statistics",NUM_EID_STATS,0.5,NUM_EID_STATS+0.5);
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
		_heid->GetXaxis()->SetBinLabel(EID_ECIN_MIN,"ec_ei > ECmin");
		_heid->GetXaxis()->SetBinLabel(EID_EC_FID,"pass ECfid");
		_heid->GetXaxis()->SetBinLabel(EID_ZVTX,"pass zvrtx");
		_heid->GetXaxis()->SetBinLabel(EID_SF,"pass sf");
		//! EFID
		_hefid=new TH1D("hefid","EFID statistics",NUM_EFID_STATS,0.5,NUM_EFID_STATS+0.5);
		_hefid->GetXaxis()->SetBinLabel(EFID_TOT,"total");
		_hefid->GetXaxis()->SetBinLabel(EFID_IN,"infid");
		//! PID
		_hpid=new TH1D("hpid","PID statistics",NUM_PID_STATS,0.5,NUM_PID_STATS+0.5);
		_hpid->GetXaxis()->SetBinLabel(PID_TOT,"total");
		_hpid->GetXaxis()->SetBinLabel(PID_PASS,"pass");
		//! pcorr
		if (_dtyp=="exp"){
			//!pcorr
			_hpcorr_dpVp=new TH2D("hdpVp","#Deltap vs. p",550,0,5.5,160,-0.08,0.08);
			_hpcorr_dcx=new TH1D("hdcx", "#Deltacx",60,-0.01,0.01);
			_hpcorr_dcy=new TH1D("hdcy", "#Deltacy",60,-0.01,0.01);
			_hpcorr_dcz=new TH1D("hdcz", "#Deltacz",60,-0.01,0.01);
			_hpcorr_dp=new TH1D("hdp", "#Deltap", 160,-0.08,0.08);
		}
	}
	//! delast
	_hW=new TH1D("hW","hW",250,0,5);
	_helast=new TH1D("helast","helast",100,0.7,1.2);
	_hf=new TH1D*[6];
	_hc=new TH1D*[6];
	for (int isctr=0;isctr<6;isctr++){
		TString name;
		TString title;
		//! _hf
		name=TString::Format("hf_s%d",isctr+1);
		title=TString::Format("#phi=[%d,%d]",PHI_FULL_SECTOR[isctr][0],PHI_FULL_SECTOR[isctr][1]);
		_hf[isctr]=new TH1D(name,title,_NBINS,_THETA_MIN,_THETA_MAX);
		//! _hc
		name=TString::Format("hc_s%d",isctr+1);
		title=TString::Format("#phi=[%d,%d]",PHI_CENTRAL_SECTOR[isctr][0],PHI_CENTRAL_SECTOR[isctr][1]);
		_hc[isctr]=new TH1D(name,title,_NBINS,_THETA_MIN,_THETA_MAX);
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
	Info("","P_MIN_ECTH=%.3f",P_MIN_ECTH);
	for (int isctr=0;isctr<6;isctr++){
		Info("","ECIN_MIN sctr%d=%.3f",isctr+1,_ECmin[isctr]);
	}
	Info("","UMIN=%.0f,UMAX=%.0f,VMIN=%.0f,VMAX=%.0f,WMIN=%.0f,WMAX=%.0f",
		      UMIN,UMAX,VMIN,VMAX,WMIN,WMAX);
	for (int isctr=0;isctr<6;isctr++){
		Info("","_z_vtx_min_sctr%d=%.2f",isctr+1,_z_vtx_min[isctr]);
		Info("","_z_vtx_max_sctr%d=%.2f",isctr+1,_z_vtx_max[isctr]);
	}
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

	delete _pid_tool;
	delete _pcorr;

	delete[] _sf_mean;
	delete[] _sf_min;
	delete[] _sf_max;
	delete _ECmin;
	delete _z_vtx_min;
	delete _z_vtx_max;
	
	delete _hevt;
	delete _heid;
	delete _hefid;
	delete _hpid;
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
   Info("h10looper_e1f::~h10looper_e1f","Done.");
}

/*
+ The following function sets up eid cut parameters. 
+ The cuts that can be directly set up, for example the ones that are the same
for exp and sim and are a simple number are directly set up constant.h.
+ Others that are different for exp and sim and require complex Objects, like the SF cut, are set up here.
*/
void h10looper_e1f::setup_eid_cutpars(TString dtyp)
{
	//! minimum momentum that satisfies EC threshold
	//! directly set in constants.h

	//! ECin cut pars
	//! pars obtained from constants.h
	_ECmin=new Float_t[6];
	for (int isctr=0;isctr<6;isctr++){
		if      (dtyp=="exp") _ECmin[isctr]=ECIN_MIN_EXP[isctr];
		else if (dtyp=="sim") _ECmin[isctr]=ECIN_MIN_SIM[isctr];
	}

	//! EC fid cut pars (as per MG and EP, and therefore as per "E1F run group"?)
	//! directly setup in constant.h

	//! z-vertex cut
	_z_vtx_min=new Float_t[6];
	_z_vtx_max=new Float_t[6];
	for (int isctr=0;isctr<6;isctr++){
		if(dtyp=="exp"){
			_z_vtx_min[isctr]=ZVTX_MIN_EXP[isctr];
			_z_vtx_max[isctr]=ZVTX_MAX_EXP[isctr];
		}else if (dtyp=="sim"){
			_z_vtx_min[isctr]=ZVTX_MIN_SIM[isctr];
			_z_vtx_max[isctr]=ZVTX_MAX_SIM[isctr];
		} 
	}

	//! SF cut
	//! pars obtained from constants.h
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
					if ( !_use_proton || (_use_proton && found_proton()) ){
						if (_W>_W_CUT_MIN && _W<_W_CUT_MAX){
							_hevt->Fill(EVT_ELSTC);
							make_delast();
						}
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

	//! prepare input for lvl_1 cut
	int gprt=gpart;
	int chrg=q[0];
	bool hitDC=dc[0]>0,hitCC=cc[0]>0,hitSC=sc[0]>0,hitEC=ec[0]>0;
	int idxDC=dc[0]-1,idxCC=cc[0]-1,idxSC=sc[0]-1,idxEC=ec[0]-1;
	int stt=stat[0], dc_stt=dc_stat[idxDC];
	
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
					if ( (!_use_SChit) || (_use_SChit && hitSC) ){
						_heid->Fill(EID_HIT_SC);
						if (hitEC){
							_heid->Fill(EID_HIT_EC);
							if (stt>0){
								_heid->Fill(EID_STAT);
								if ( (!_use_dc_stat) ||(_use_dc_stat && dc_stt>0) ){
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
		if (pass_p_min_ECth()){
			_heid->Fill(EID_P_MIN_ECTH);
			if ( (!_use_cut_ECin_min) || (_use_cut_ECin_min && pass_ECin_min()) ){
				_heid->Fill(EID_ECIN_MIN);
				if ( (!_use_cut_ECfid) || (_use_cut_ECfid && pass_ECfid()) ){
					_heid->Fill(EID_EC_FID);
					if( (!_use_cut_zvtx) || (_use_cut_zvtx && pass_zvtx()) ){
						_heid->Fill(EID_ZVTX);
						if (pass_sf()){
							_heid->Fill(EID_SF);
							ret=kTRUE;
						} 
					}
				}
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
	bool ret=kFALSE;
	if (_use_ep_efid){
		ret=inFid(mom,_theta,_phi,sector);
	}else{
		ret=Cuts::Fiducial(id, mom, _theta, _phi, sector, paddle);
	}
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

/*
found_proton() alogrithm assumes that there is only 
1 proton in an event. Therefore once the requirements
for a track being a proton is fulfilled, the function returns.
*/
bool h10looper_e1f::found_proton(){
	bool ret=kFALSE;
	_hpid->Fill(PID_TOT);
	//! Directly measured electron quantities
	float l_e=sc_r[sc[0]-1];
	float t_e=sc_t[sc[0]-1];
	float t_off=t_e-(l_e/SOL);
	//! Now loop over gpart to find proton
	for (int i=1;i<gpart;i++){
		int chrg=q[i];
		bool hitDC=dc[i]>0,hitSC=sc[i]>0;
		int stt=stat[i], dc_stt=dc_stat[dc[i]-1];

		if(chrg==1 && hitDC && hitSC){
			//! Get time of flight for track
			float t=sc_t[sc[i]-1];

			//! Now calculate dt using :
			//! + length of track and 
			//! + beta of track under assumption that it is a proton
			//! that particle is a proton
			float l=sc_r[sc[i]-1];
			//! Calculate beta under assumption that particle is proton
			float mom=p[i];
			float b_p=TMath::Sqrt((mom*mom)/(MASS_P*MASS_P+mom*mom));
			//! Now calculate dt
			float dt=l/(b_p*SOL)+t_off- t;

			//! Using dt and p, identify if particle is a proton
			if (_pid_tool->is_proton(dt,mom)){
				ret=kTRUE;
				_hpid->Fill(PID_PASS);
				break;
			}
		}
	}
	return ret;
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

/*****************************************************
[07-09-15]
The following function, taken directly as sent to my by Evan,
takes the global x,y,z coordinates of the EC (ech_x,ech_y,ech_z)
and fills calculates the local U,V,W coordinates of the EC.
*****************************************************/
void h10looper_e1f::GetUVW(float xyz[3], float uvw[3]) {
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
  phi = TMath::ATan2(xyz[Y], xyz[X]) * 57.29578;
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

bool h10looper_e1f::pass_sf()
{
	//! Get info to make cut
	float mom=p[0];
	int idxEC=ec[0]-1;
	int sctr=ec_sect[idxEC];
	//! correct etot. etot=max(etot,ec_ei_ec_eo)
	if (_use_corr_sf_etot){
		etot[idxEC]=etot[idxEC] > ec_ei[idxEC]+ec_eo[idxEC] ? etot[idxEC] : ec_ei[idxEC]+ec_eo[idxEC];
	}
	float sf=etot[idxEC]/mom;

	//! make cut
	float sf_min=_sf_min[sctr-1]->Eval(mom);
	float sf_max=_sf_max[sctr-1]->Eval(mom);
	return ( sf>sf_min && sf<sf_max );
}

bool h10looper_e1f::pass_ECfid(){
	int idxEC=ec[0]-1;
	Float_t xyz[3]={ech_x[idxEC],ech_y[idxEC],ech_z[idxEC]};
	Float_t uvw[3]={0,0,0};
	GetUVW(xyz,uvw);
	enum { U, V, W };
	bool passU= uvw[U]>UMIN && uvw[U]<UMAX;
	bool passV= uvw[V]>VMIN && uvw[V]<VMAX;
	bool passW= uvw[W]>WMIN && uvw[W]<WMAX;
	return passU && passV && passW;
}

bool h10looper_e1f::pass_p_min_ECth(){
	return (p[0]>P_MIN_ECTH);
}

bool h10looper_e1f::pass_ECin_min(){
	int idxEC=ec[0]-1;
	int sctr=ec_sect[idxEC];
	return (ec_ei[idxEC]>_ECmin[sctr-1]);
}

bool h10looper_e1f::pass_zvtx(){
	float zvtx=vz[0];
	int idxEC=ec[0]-1;
	int sctr=ec_sect[idxEC];
	
	return (zvtx>_z_vtx_min[sctr-1] && zvtx<_z_vtx_max[sctr-1]);
}

void h10looper_e1f::setup_adtnl_cut_opts(TString adtnl_cut_opt){
	//! default values
	//! cuts to be made in addition to 'dflt'
    //! eid: new cuts and corrections
    _use_cut_ECin_min=kFALSE;
	_use_cut_ECfid=kFALSE;
	_use_cut_zvtx=kFALSE;
	_use_corr_sf_etot=kFALSE;
	//! eid: to match Isupov
	_use_SChit=kTRUE;
	_use_dc_stat=kTRUE;
	//! Evans's EFID
	_use_ep_efid=kFALSE;
	//! proton
	_use_proton=kFALSE;

	//! Now set as per adtnl_cut_opt
	if (adtnl_cut_opt.Contains("1:")) _use_cut_ECin_min=kTRUE;
	if (adtnl_cut_opt.Contains("2:")) _use_cut_ECfid=kTRUE;
	if (adtnl_cut_opt.Contains("3:")) _use_cut_zvtx=kTRUE;
	if (adtnl_cut_opt.Contains("4:")) _use_corr_sf_etot=kTRUE;
	if (adtnl_cut_opt.Contains("5:")) _use_SChit=kFALSE;
	if (adtnl_cut_opt.Contains("6:")) _use_dc_stat=kFALSE;
	if (adtnl_cut_opt.Contains("7:")) _use_ep_efid=kTRUE;
	if (adtnl_cut_opt.Contains("8:")) _use_proton=kTRUE;

	Info("h10looper_e1f::setup_adtnl_cut_opts","***adtnl_cut_opt=%s***",adtnl_cut_opt.Data());
	Info("h10looper_e1f::setup_adtnl_cut_opts","The following cuts-corrections will be made in addition to \'dflt\':");
	//! eid: new cuts and corrections
	if (_use_cut_ECin_min) Info("","use_cut_ECin_min");
	if (_use_cut_ECfid) Info("","use_cut_ECfid");
	if (_use_cut_zvtx) Info("","use_cut_zvtx");
	if (_use_corr_sf_etot) Info("","use_corr_sf_etot");
	//! eid: to match Isupov
	if (!_use_SChit) Info("","not_use_SChit");
	if (!_use_dc_stat) Info("","not_use_dc_stat");
	//! Evans's EFID
	if(_use_ep_efid) Info("","use_ep_efid");
	//! proton 
	if (_use_proton) Info("","use_proton");
	Info("","***********");
}
