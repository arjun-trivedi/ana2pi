#define h10looper_2pi_cxx
#include "h10looper_2pi.h"

#include "h8_bng.h"

#include "fidfuncs_hadrons.C" //! EP_PFID
#include "wrpr_cut_fid_e16.h" //! EI_PFID

#include <TLorentzRotation.h>

h10looper_2pi::h10looper_2pi(TString h10type, TChain* h10chain,
	                         TString cutsncors,
	                         TString fout_name, Long64_t nentries,
						     TString adtnl_opts) : h10looper_e1f(h10type,h10chain,cutsncors,fout_name,nentries,adtnl_opts)
{
	Info("h10looper_2pi::h10looper_2pi","Setting up h10looper_2pi...\n");

	//! output objects
	//! pfid
	if (_do_pfid){
		_fout->mkdir("pfid")->cd();
		_hpfid=new TH1D("hpfid","PFID statistics",NUM_PFID_STATS,0.5,NUM_PFID_STATS+0.5);
		_hpfid->GetXaxis()->SetBinLabel(PFID_TOT,"total");
		_hpfid->GetXaxis()->SetBinLabel(PFID_P_IN,"p infid");
		_hpfid->GetXaxis()->SetBinLabel(PFID_PIP_IN,"#pi^{+} infid");
		_hpfid->GetXaxis()->SetBinLabel(PFID_P_AND_PIP_IN,"p + #pi^{+} infid");
		_hpfid->SetMinimum(0);
		//! phi vs. theta hists
		_hpfid_p=new TH2F*[2];
		_hpfid_pip=new TH2F*[2];
		for (int i=0;i<2;i++){
			TString name_sfx;
			if      (i==0) name_sfx="prec";
			else if (i==1) name_sfx="pstc";
			_hpfid_p[i]  = new TH2F(TString::Format("h_p_phiVtheta_%s",name_sfx.Data()),  "#phi vs. #theta for p",  100,0,60,  100,-30,330);
			_hpfid_pip[i]= new TH2F(TString::Format("h_pip_phiVtheta_%s",name_sfx.Data()),"#phi vs. #theta for pip",100,0,120, 100,-30,330);
		}
	}
	//!eff
	if (_do_eff){
		_fout->mkdir("eff")->cd();
		_heff=new TH1D("heff","EFF statistics",NUM_EFF_STATS,0.5,NUM_EFF_STATS+0.5);
		_heff->GetXaxis()->SetBinLabel(EFF_TOT,"total");
		_heff->GetXaxis()->SetBinLabel(EFF_E_PASS,"e^{-} ineff");
		_heff->GetXaxis()->SetBinLabel(EFF_P_PASS,"p ineff");
		_heff->GetXaxis()->SetBinLabel(EFF_PIP_PASS,"#pi^{+} ineff");
		_heff->GetXaxis()->SetBinLabel(EFF_E_AND_P_AND_PIP_PASS,"e+p+#pi^{+} ineff");
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
	//!scpd
	if (_do_scpd){
		_fout->mkdir("scpd")->cd();
		_hscpd=new TH1D("hscpd","SCPD statistics",NUM_SCPD_STATS,0.5,NUM_SCPD_STATS+0.5);
		_hscpd->GetXaxis()->SetBinLabel(SCPD_TOT,"total");
		_hscpd->GetXaxis()->SetBinLabel(SCPD_E_PASS,"e^{-} inscpd");
		_hscpd->GetXaxis()->SetBinLabel(SCPD_P_PASS,"p inscpd");
		_hscpd->GetXaxis()->SetBinLabel(SCPD_PIP_PASS,"#pi^{+} inscpd");
		_hscpd->GetXaxis()->SetBinLabel(SCPD_E_AND_P_AND_PIP_PASS,"e+p+#pi^{+} inscpd");
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
					if      (_eff_prtcl_names[i]=="e")   _hthetavp2[i][j][k]=new TH2D(name,name,100,0,5,100,0,60);
					else if (_eff_prtcl_names[i]=="p")   _hthetavp2[i][j][k]=new TH2D(name,name,100,0,4,100,0,60);
					else if (_eff_prtcl_names[i]=="pip") _hthetavp2[i][j][k]=new TH2D(name,name,100,0,3,100,0,120);	
				}
			}
		}
	}
	//! d2pi
	if (_do_evtsel_2pi){
		setup_d2pi();
	}
	if(_make_h10_skim_SS){
		//! + h10 is created directly under 'fout:/' for technical reasons
		//!   (since the program expects to find the h10 tree in the root dir.)
		//! + Directory 'copyh10' is created only for book-keeping reasons
		_fout->mkdir("copyh10");
		_fout->cd();
		_th10copy = (TTree*)fChain->GetTree()->CloneTree(0);
		set_h10_SEB_BranchStatus(_th10copy);

		/* 
		[02-28-16] hack-pcorr for '_make_h10_skim_SS'
		+ 'pcorr' is not a part of proc-chain when making h10_skim_SS: its usage is reserved
		   for when the h10 variables are to be updated with the corrected momenta which are
		   then used to updated electron and hadron kinematics via 'set_ekin()' and 'set_hkin()'
		+ However since in 'evtsel', that is used in the proc-chain, the corrected momenta have
		  to be used but h10 is not updated, I had to hack 'pcorr' in a manner defined under the header 
		  "hack-pcorr for '_make_h10_skim_SS'" in h10looper_e1f(2pi).h
		+ In this hack-pcorr method, 'pcorr' is put inside the 'd2pi' folder  
		*/

		if (_dtyp=="exp"){
			_pcorr_prtcl_names=new TString[3];
			_pcorr_prtcl_names[0]="e";
			_pcorr_prtcl_names[1]="p";
			_pcorr_prtcl_names[2]="pip";

			TDirectory* dir_d2pi=_fout->GetDirectory("d2pi");
			dir_d2pi->mkdir("pcorr")->cd();
			_hpcorr_dpVp=new TH2D*[3];
			_hpcorr_dcx=new TH1D*[3];
			_hpcorr_dcy=new TH1D*[3];
			_hpcorr_dcz=new TH1D*[3];
			_hpcorr_dp=new TH1D*[3];
			for (int i=0;i<3;i++){
				TString sfx=_pcorr_prtcl_names[i].Data();
				_hpcorr_dpVp[i]=new TH2D(TString::Format("hdpVp_%s",sfx.Data()),"#Deltap vs. p",550,0,5.5,160,-0.08,0.08);
				_hpcorr_dcx[i]=new TH1D(TString::Format("hdcx_%s",sfx.Data()), "#Deltacx",60,-0.01,0.01);
				_hpcorr_dcy[i]=new TH1D(TString::Format("hdcy_%s",sfx.Data()), "#Deltacy",60,-0.01,0.01);
				_hpcorr_dcz[i]=new TH1D(TString::Format("hdcz_%s",sfx.Data()), "#Deltacz",60,-0.01,0.01);
				_hpcorr_dp[i]=new TH1D(TString::Format("hdp_%s",sfx.Data()), "#Deltap", 160,-0.08,0.08);
			}
		}
	
	}


	//! Hadron kinematics
	//! Lab Frame
	_lvP.SetXYZT(0,0,0,0);
	_lvPip.SetXYZT(0,0,0,0);
	_lvPim.SetXYZT(0,0,0,0);
	//! Varset
	_lvQCMS.SetXYZT(0,0,0,0);
    _lvP0CMS.SetXYZT(0,0,0,0);
    _lvPCMS.SetXYZT(0,0,0,0);
    _lvPipCMS.SetXYZT(0,0,0,0);
    _lvPimCMS.SetXYZT(0,0,0,0);
	//! Lab frame
	_p_p=  _theta_p=  _phi_p=0;
    _p_pip=_theta_pip=_phi_pip=0;
    _p_pim=_theta_pim=_phi_pim=0;
    //! Varsets
    _M_ppip=_M_ppim=_M_pippim=0;
    _p_cms_p=  _theta_cms_p=  _phi_cms_p=0;
    _p_cms_pip=_theta_cms_pip=_phi_cms_pip=0;
    _p_cms_pim=_theta_cms_pim=_phi_cms_pim=0;
    _alpha_1=_alpha_2=_alpha_3=0;
	
	Info("h10looper_2pi::h10looper_2pi","Following information was setup:");
	Info("","---------------------------------------");
	if (_do_pfid){
		Info("","*** pfid pars ***");
		Info("","Printing not implemented here!");
	}
	if (_do_evtsel_2pi){
		Info("","*** evtsel_2pi pars ***");
		Info("","_mm2ppip_l=%f,_mm2ppip_h=%f",_mm2ppip_l,_mm2ppip_h);
	}
	Info("","---------------------------------------");
	Info("h10looper_2pi::h10looper_2pi","Done setting up h10looper_2pi\n");
	
}

h10looper_2pi::~h10looper_2pi()
{
	Info("h10looper_2pi::~h10looper_2pi","");
	
	delete _hpfid;
	delete[] _hpfid_p;
	delete[] _hpfid_pip;
	delete _heff;
	delete[] _eff_prtcl_names;
	delete[] _hthetavp;
	delete _hscpd;
	delete[] _scpd_prtcl_names;
	delete[] _hpdl;
	delete[] _hthetavp2;

	delete _hq2w_prec;
	delete _hq2w_pstc;
	delete _hmm2_prec_fW;
	delete _hmm_prec_fW;
	delete _hmm2_pstc_fW;
	delete _hmm_pstc_fW;
	delete[] _hmm2_prec;
	delete[] _hmm_prec;
	delete[] _hmm2_pstc;
	delete[] _hmm_pstc;
	delete[] _h8;

	//![02-28-16] After adding hack-pcorr
	if (_make_h10_skim_SS){
		delete[] _pcorr_prtcl_names;
		delete[] _hpcorr_dpVp;
		delete[] _hpcorr_dcx;
   		delete[] _hpcorr_dcy;
   		delete[] _hpcorr_dcz;
   		delete[] _hpcorr_dp;
    }

	Info("h10looper_2pi::~h10looper_2pi","Done.");
}

void h10looper_2pi::Loop(){
	Info("h10looper_2pi::Loop()","");

	if (fChain == 0) return;

	//!set_h10_SEB_BranchStatus()
	if (_seq=="recon"){
		set_h10_SEB_BranchStatus(fChain);
	}

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<_nentries_to_proc;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		if ( (jentry+1)%100000==0 ) {
			Info("h10looper_2pi::Loop()","Processing event %llu",jentry+1);
			Info("h10looper_2pi::Loop()","Fraction of events processed %.2f%%",((float)(jentry+1)/_nentries_to_proc)*100);
		}

		//! Reconcile h10 Branch binding for e16:exp
		//! or
		//! if _do_reconcile specified by user (for example for sim-e16-EI)
		if (_do_reconcile || (_expt=="e16" && _dtyp=="exp")) {
			Reconcile();
		}

		//! Begin processing event
		reset_ekin();
		reset_hkin();
		if (_seq=="recon"){
			_hevt->Fill(EVT_TRG);
			//! EID
			if (_do_eid)_heid->Fill(EID_TRG);
			//std::cout<<"eid"<<std::endl;
			if ( !_do_eid  || (_do_eid && pass_eid()) ){
				if (_do_eid)_heid->Fill(EID_E);
				_hevt->Fill(EVT_E);
				set_ekin();
				//! EFID
				if (_do_efid) _hefid->Fill(EFID_TOT);
				//std::cout<<"efid"<<std::endl;
				if ( !_do_efid || (_do_efid && pass_efid()) ){
					if (_do_efid) _hefid->Fill(EFID_IN);
					_hevt->Fill(EVT_E_INFID);

					//! make_h10_skim_e
					if (_make_h10_skim_e){
						//std::cout<<"h10-skim-e enter"<<std::endl;
						_th10copy->Fill();
						//std::cout<<"h10-skim-e leave"<<std::endl;
						continue;
					}

					/*//! pcorr
					if (_do_pcorr){
						//mom_corr_electron();
						do_pcorr();
						//std::cout<<"pcorr"<<std::endl;
						set_ekin();
					}*/
					//! PID (top2')
					if (_do_pid){
						_hpid->Fill(PID_TOT);
						fill_pid_prec_hists();
					}
					//std::cout<<"pid"<<std::endl;
					int h10idx_p=0,h10idx_pip=0;
					if ( !_do_pid || (_do_pid && pass_pid(h10idx_p,h10idx_pip)) ){
						if(_do_pid) {
							_hpid->Fill(PID_P_AND_PIP_FOUND);
							fill_pid_pstc_hists(h10idx_p, h10idx_pip);
						}
						_hevt->Fill(EVT_P_PIP);
						set_hkin(h10idx_p, h10idx_pip);
						//! PFID (top2')
						if (_do_pfid) _hpfid->Fill(PFID_TOT);
						//std::cout<<"pfid"<<std::endl;
						if ( !_do_pfid || (_do_pfid && proton_infid() && pip_infid()) ){
							if (_do_pfid) _hpfid->Fill(PFID_P_AND_PIP_IN);
							_hevt->Fill(EVT_P_PIP_INFID);

							//! pcorr
							if (_do_pcorr){
								//mom_corr_electron();
								do_pcorr();
								//std::cout<<"pcorr"<<std::endl;
								set_ekin();
								set_hkin(h10idx_p, h10idx_pip);
							}
														
							//! EFF (top2')
							if (_do_eff) _heff->Fill(EFF_TOT);
							if (!_do_eff || (_do_eff && pass_eff()) ){
								if (_do_eff) _heff->Fill(EFF_E_AND_P_AND_PIP_PASS);
								_hevt->Fill(EVT_INEFF);

								//!SCPD
								if (_do_scpd) _hscpd->Fill(SCPD_TOT);
								if (!_do_scpd || (_do_scpd && pass_scpd()) ){
									if (_do_scpd) _hscpd->Fill(SCPD_E_AND_P_AND_PIP_PASS);
									_hevt->Fill(EVT_INSCPD);

									/*//! pcorr
									if (_do_pcorr){
										//mom_corr_electron();
										do_pcorr();
										//std::cout<<"pcorr"<<std::endl;
										set_ekin();
										set_hkin(h10idx_p, h10idx_pip);
									}*/
								
									//! Event selection
									if (_do_evtsel_2pi){
										//! Q2-W kinematic cut
										_hq2w_prec->Fill(_W,_Q2);
										//! hack-pcorr
										if (_make_h10_skim_SS && _dtyp=="exp"){//! then do hack-pcorr
											TLorentzVector lv_corr_E,lv_corr_P,lv_corr_Pip;
											lv_corr_E.SetXYZT(0,0,0,0);
											lv_corr_P.SetXYZT(0,0,0,0);
											lv_corr_Pip.SetXYZT(0,0,0,0);
											//std::cout<<"here:start"<<std::endl;
											do_pcorr_but_not_update_h10(lv_corr_E,lv_corr_P,lv_corr_Pip);
											set_ekin_use_passed_lv(lv_corr_E);
											set_hkin_use_passed_lv(lv_corr_P,lv_corr_Pip,h10idx_p,h10idx_pip);
											//std::cout<<"here:done"<<std::endl;
										}
										//! [11-23-15] see h8_bng.h for details if (!(_Q2>=1.25 && _Q2<5.25 && _W>1.300 && _W<2.125)) continue;
										if (!(_Q2>=_Q2_MIN && _Q2<_Q2_MAX && _W>_W_MIN && _W<_W_MAX)) continue;
										_hq2w_pstc->Fill(_W,_Q2);
										_hevt->Fill(EVT_Q2W_KIN_PASS);

										//! Topology selection
										//! + Till this point, pid is done for t2'
										//! + If _use_t2 is True, then t2' -> t2 i.e. exclusive p,pip,pim_missing
										if (_use_t2 && found_hadron("pim")>=1) continue;
			
										//! MM cut
										//! + NOTE that _lvPim is recons'd with no knowledge of pim!
										//!    + _lvPim=_lvW-(_lvP+_lvPip)
										//! + Therefore, its Mag2() and Mag() can be treated as MM2 and MM
										float mm2ppip=_lvPim.Mag2();
										float mmppip=_lvPim.Mag();
										_hmm2_prec_fW->Fill(mm2ppip);
										_hmm_prec_fW->Fill(mmppip);
										int iw=GetCrsWBinIdx(_W);
										_hmm2_prec[iw]->Fill(mm2ppip);
										_hmm_prec[iw]->Fill(mmppip);
										if (mm2ppip>_mm2ppip_l && mm2ppip<_mm2ppip_h){
											_hevt->Fill(EVT_2PI);
											_hmm2_pstc_fW->Fill(mm2ppip);
											_hmm_pstc_fW->Fill(mmppip);
											_hmm2_pstc[iw]->Fill(mm2ppip);
											_hmm_pstc[iw]->Fill(mmppip);

											//! make_h10_skim_SS
											if (_make_h10_skim_SS){
												//std::cout<<"h10-skim-SS enter"<<std::endl;
												_th10copy->Fill();
												//std::cout<<"h10-skim-SS leave"<<std::endl;
												continue;
											}

											fill_h8();
										}//! MM2-cut
									}//! EVTSEL_2PI
								}//!SCPD
							}//!EFF
						}//!PFID
					}//!PID
				}//!EFID
			}//!EID
		}else if (_seq=="thrown"){
			_hevt->Fill(EVT_TRG);
			set_ekin();
			set_hkin();
			
			//! Event selection

			//! Q2-W kinematic cut
			_hq2w_prec->Fill(_W,_Q2);
			if (!(_Q2>=_Q2_MIN && _Q2<_Q2_MAX && _W>_W_MIN && _W<_W_MAX)) continue;
			_hevt->Fill(EVT_Q2W_KIN_PASS);
			_hq2w_pstc->Fill(_W,_Q2);

			_hevt->Fill(EVT_2PI);
			fill_h8();
		}
	}//end event loop

	//! Write _fout
	_fout->Write();
}//end Loop()

void h10looper_2pi::reset_hkin(){
	//! Lab Frame
	_lvP.SetXYZT(0,0,0,0);
	_lvPip.SetXYZT(0,0,0,0);
	_lvPim.SetXYZT(0,0,0,0);
	//! Varset
	_lvQCMS.SetXYZT(0,0,0,0);
    _lvP0CMS.SetXYZT(0,0,0,0);
    _lvPCMS.SetXYZT(0,0,0,0);
    _lvPipCMS.SetXYZT(0,0,0,0);
    _lvPimCMS.SetXYZT(0,0,0,0);
	//! Lab frame
	_p_p=  _theta_p=  _phi_p=0;
    _p_pip=_theta_pip=_phi_pip=0;
    _p_pim=_theta_pim=_phi_pim=0;
    //! sector and sc_pd
    _sector_p=_sc_pd_p=0;
    _sector_pip=_sc_pd_pip=0;
    _sector_pim=_sc_pd_pim=0;
    //! Varset kinematics
    _M_ppip=_M_ppim=_M_pippim=0;
    _p_cms_p=  _theta_cms_p=  _phi_cms_p=0;
    _p_cms_pip=_theta_cms_pip=_phi_cms_pip=0;
    _p_cms_pim=_theta_cms_pim=_phi_cms_pim=0;
    _alpha_1=_alpha_2=_alpha_3=0;
    //! h10idx
    _h10idx_p=_h10idx_pip=_h10idx_pim=-1;
}

/*
set_hkin() implemented as per top2' logic
*/
void h10looper_2pi::set_hkin(int h10idx_p/*=-1*/, int h10idx_pip/*=-1*/){
	if (_seq=="recon"){
		//! First set _lvP
		float mom_p = p[h10idx_p];
		float px_p = mom_p*cx[h10idx_p];
		float py_p = mom_p*cy[h10idx_p];
		float pz_p = mom_p*cz[h10idx_p];
		float e_p = TMath::Sqrt(mom_p*mom_p+MASS_P*MASS_P);
		_lvP.SetPxPyPzE(px_p,py_p,pz_p,e_p);

		//! Now set _lvPip
		float mom_pip = p[h10idx_pip];
		float px_pip = mom_pip*cx[h10idx_pip];
		float py_pip = mom_pip*cy[h10idx_pip];
		float pz_pip = mom_pip*cz[h10idx_pip];
		float e_pip = TMath::Sqrt(mom_pip*mom_pip+MASS_PIP*MASS_PIP);
		_lvPip.SetPxPyPzE(px_pip,py_pip,pz_pip,e_pip);

		//! Now set _lvPim
		_lvPim = (_lvW-(_lvP+_lvPip));
	}else if (_seq=="thrown"){
		for (int idx=0;idx<mcnentr;idx++) {
			int id=mcid[idx];
			float mom=mcp[idx];
			float theta=mctheta[idx]*TMath::DegToRad();
			float phi=mcphi[idx]*TMath::DegToRad();
			float mass=mcm[idx];
			float energy=TMath::Sqrt(mom*mom+mass*mass);
			float pz=mom*TMath::Cos(theta);
			float py=mom*TMath::Sin(phi)*TMath::Sin(theta);
			float px=mom*TMath::Cos(phi)*TMath::Sin(theta);
			switch(id) {
				case PROTON:
					_lvP.SetPxPyPzE(px,py,pz,energy);
					break;
				case PIP:
					_lvPip.SetPxPyPzE(px,py,pz,energy);
					break;
				case PIM:
					_lvPim.SetPxPyPzE(px,py,pz,energy);
					break;
				default:
					break;
			}
		}
	}
	
	//! Set up Lab frame components of kin.
	_p_p=_lvP.P();
	_p_pip=_lvPip.P();
	_p_pim=_lvPim.P();
	_theta_p=_lvP.Theta()*TMath::RadToDeg();
	_theta_pip=_lvPip.Theta()*TMath::RadToDeg();
	_theta_pim=_lvPim.Theta()*TMath::RadToDeg();
	//! phi: 1st get in range [-180,180]
	float phi_p=_lvP.Phi()*TMath::RadToDeg();// [-180,180]
	float phi_pip=_lvPip.Phi()*TMath::RadToDeg();// [-180,180]
	float phi_pim=_lvPim.Phi()*TMath::RadToDeg();// [-180,180]
	//! phi: now convert to [-30,330]
	_phi_p=phi_p<-30?phi_p+360:phi_p; // [-30,330]
	_phi_pip=phi_pip<-30?phi_pip+360:phi_pip; // [-30,330]
	_phi_pim=phi_pim<-30?phi_pim+360:phi_pim; // [-30,330]
	//! sector
	_sector_p=get_sector(_phi_p);
	_sector_pip=get_sector(_phi_pip);
	_sector_pim=get_sector(_phi_pim);
	//! paddle
	if (_seq=="recon"){
		_sc_pd_p=100*dc_sect[dc[h10idx_p]-1]+sc_pd[sc[h10idx_p]-1];
		_sc_pd_pip=100*dc_sect[dc[h10idx_pip]-1]+sc_pd[sc[h10idx_pip]-1];
		//!_sc_pd_pim= thus far, not identifying pim
	}

	//! + The following code taken directly from ProcD2pi::UpdateD2pi
	//!   and therein the part to calculate Varset kinematics
	//! + The only changes relate to certain variable names that have 
	//!   changed and code readability etc.

	//! Varsets kinematics
	//! Calculate rotation: taken from Evan's phys-ana-omega on 08-05-13
	TVector3 uz = _lvQ.Vect().Unit();
    TVector3 ux = (_lvE0.Vect().Cross(_lvE1.Vect())).Unit();
    ux.Rotate(-TMath::Pi()/2,uz);
    TRotation r3;// = new TRotation();
    r3.SetZAxis(uz,ux).Invert();
    //! _w and _q are in z-direction
    TVector3 boost(-1*_lvW.BoostVector());
    TLorentzRotation r4(r3); //*_boost);
    r4 *= boost; //*_3rot;
	
	_lvQCMS   = _lvQ;
	_lvP0CMS  = _lvP0;
	_lvPCMS   = _lvP;
	_lvPipCMS = _lvPip;
	_lvPimCMS = _lvPim;
	//! Rotate and Boost
	_lvQCMS.Transform(r4);
	_lvP0CMS.Transform(r4);
	_lvPCMS.Transform(r4);
	_lvPipCMS.Transform(r4);
	_lvPimCMS.Transform(r4);
	
	_M_ppip = (_lvPCMS + _lvPipCMS).Mag();
	_M_ppim = (_lvPCMS + _lvPimCMS).Mag();
	_M_pippim = (_lvPipCMS + _lvPimCMS).Mag();
	
	//! Should be directly be able to get the Theta() angle since
	//! virtual photon defines the coordinate system
	_theta_cms_p=_lvPCMS.Theta()*TMath::RadToDeg();
	_theta_cms_pip=_lvPipCMS.Theta()*TMath::RadToDeg();
	_theta_cms_pim=_lvPimCMS.Theta()*TMath::RadToDeg();
	
	_phi_cms_p=getPhi(_lvPCMS);
	_phi_cms_pip=getPhi(_lvPipCMS);
	_phi_cms_pim=getPhi(_lvPimCMS);

	//! alpha angle
	//! Following vectors are used in the calculation of alpha angle; see doc. for getAlpha()
	TVector3 G_f(0,0,0);
	TVector3 G_p(0,0,0);
	TVector3 B_f(0,0,0);
	TVector3 B_p(0,0,0);


	G_f=uz*-1; //! Gamma vector used in calculation of alpha always "follows" -uz

	//! alpha[p',pip][p,pim]
	G_p=_lvPimCMS.Vect().Unit();
	B_f=_lvPipCMS.Vect().Unit();
	B_p=G_p;
	_alpha_1=getAlpha(G_f,G_p,B_f,B_p);
	//! alpha[pip,pim][p,p']
	G_p=_lvPCMS.Vect().Unit();
	B_f=_lvPipCMS.Vect().Unit();
	B_p=G_p;
	_alpha_2=getAlpha(G_f,G_p,B_f,B_p);
	//! alpha[p',pim][p,pip]
	G_p=_lvPipCMS.Vect().Unit();
	B_f=_lvPimCMS.Vect().Unit();//_lvPCMS.Vect().Unit();
	B_p=G_p;
	_alpha_3=getAlpha(G_f,G_p,B_f,B_p);

	//! h10idx
	_h10idx_p=h10idx_p;
	_h10idx_pip=h10idx_pip;
	//_h10idx_pim=TBD;
}

//! PID only for top2':
//!   + only p and pip need to be detected
//!   + pim is always reconstructed using lvPim=lvW-(lvP+lvPip)
//!     i.e. it does not matter if it is detected or missing
bool h10looper_2pi::pass_pid(int& h10idx_p, int& h10idx_pip){
	bool ret=kFALSE;
	if ( !_use_gpart_pid || (_use_gpart_pid && (gpart==3 || gpart==4)) ){
		_hpid->Fill(PID_GPART_PASS);
		h10idx_p=found_hadron("p");
		h10idx_pip=found_hadron("pip");
		/*if (h10idx_p>=1) {
			_hpid->Fill(PID_P_FOUND);
		}
		if (h10idx_pip>=1) {
			_hpid->Fill(PID_PIP_FOUND);
		}*/
		if (h10idx_p>=1 && h10idx_pip>=1 && h10idx_p!=h10idx_pip){
			ret=kTRUE;
		}
	}
	return ret;
}

//! PID only for top2'
bool h10looper_2pi::pass_pfid(){
	bool ret=kFALSE;
	if (proton_infid() && pip_infid()){
		ret=kTRUE;
	}
	return ret;
}

bool h10looper_2pi::proton_infid(){
	bool ret=kFALSE;
	
	//!prec
	_hpfid_p[0]->Fill(_theta_p,_phi_p);

	int sctr_p=get_sector(_phi_p);
	if (_expt=="e1f"){
		if (_use_ep_pfid){
			//! + The last argument is momentum, but cut is independent of it
        	//! + Evan obtain cut pars for exp and sim, however, told me that exp pars have to be applied for both.
        	TF1 f_l=fPhiFid_hdrn_l_mod(PROTON,"exp",sctr_p,1);//! The last argument is mom, but cut is independent of it
        	TF1 f_h=fPhiFid_hdrn_h_mod(PROTON,"exp",sctr_p,1);//! The last argument is mom, but cut is independent of it
        	TLine l=lt0(PROTON, _dtyp);

        	float theta_min=l.GetX1();
        	float phi_min=f_l.Eval(_theta_p);
        	float phi_max=f_h.Eval(_theta_p);

        	if ( (_theta_p > theta_min) && (_phi_p > phi_min) && (_phi_p < phi_max) ){
        		ret=kTRUE;
        	}
    	}else{
			ret=Fiducial_e16_hdrn(_theta_p,_phi_p,sctr_p);
		}
	}else if (_expt=="e16"){
		ret=Fiducial_e16_hdrn(_theta_p,_phi_p,sctr_p);
	}

	if (ret==kTRUE){
		_hpfid->Fill(PFID_P_IN);	
		//!pstc
		_hpfid_p[1]->Fill(_theta_p,_phi_p);
	} 
	return ret;
}

bool h10looper_2pi::pip_infid(){
	bool ret=kFALSE;

	//!prec
	_hpfid_pip[0]->Fill(_theta_pip,_phi_pip);
	
	int sctr_pip=get_sector(_phi_pip);
	if (_expt=="e1f"){
		if (_use_ep_pfid){
			//! + The last argument is momentum, but cut is independent of it
        	//! + Evan obtain cut pars for exp and sim, however, told me that exp pars have to be applied for both.
        	TF1 f_l=fPhiFid_hdrn_l_mod(PIP,"exp",sctr_pip,1);//! The last argument is mom, but cut is independent of it
        	TF1 f_h=fPhiFid_hdrn_h_mod(PIP,"exp",sctr_pip,1);//! The last argument is mom, but cut is independent of it
        	TLine l=lt0(PIP, _dtyp);

        	float theta_min=l.GetX1();
        	float phi_min=f_l.Eval(_theta_pip);
        	float phi_max=f_h.Eval(_theta_pip);

        	if ( (_theta_pip > theta_min) && (_phi_pip > phi_min) && (_phi_pip < phi_max) ){
        		ret=kTRUE;
        	}
    	}else{
			ret=Fiducial_e16_hdrn(_theta_pip,_phi_pip,sctr_pip);
		}
	}else if (_expt=="e16"){
		ret=Fiducial_e16_hdrn(_theta_pip,_phi_pip,sctr_pip);
	}

	if (ret==kTRUE){
		 _hpfid->Fill(PFID_PIP_IN);
		 //!pstc
		_hpfid_pip[1]->Fill(_theta_pip,_phi_pip);
	}
	return ret;
}

float h10looper_2pi::getPhi(TLorentzVector lv) {
	float retVal = 0;
	retVal = TMath::RadToDeg()*invTan(lv.Py(), lv.Px());
	return retVal;
}

float h10looper_2pi::invTan(float y, float x){
	float retVal = 0; 
	if (x > 0 && y > 0)  retVal = TMath::ATan(y/x);          //1st Quad.
	if (x < 0 && y > 0)  retVal = TMath::ATan(y/x) + TMath::Pi();   //2nd Quad
	if (x < 0 && y < 0)  retVal = TMath::ATan(y/x) + TMath::Pi();   //3rd Quad
	if (x > 0 && y < 0)  retVal = TMath::ATan(y/x) + 2*TMath::Pi(); //4th Quad 
	if (x == 0 && y > 0) retVal = TMath::Pi()/2;
	if (x == 0 && y < 0) retVal = 3*TMath::Pi()/2; 
	return retVal;  
}

float h10looper_2pi::getAlpha(TVector3 uv_Gf,TVector3 uv_Gp,TVector3 uv_Bf,TVector3 uv_Bp){
	//! step 1.
	float aG=TMath::Sqrt( 1/( 1-TMath::Power(uv_Gf.Dot(uv_Gp),2) ) );
	float bG=-uv_Gf.Dot(uv_Gp)*aG;
	TVector3 v_G=aG*uv_Gf + bG*uv_Gp;
	/*Info("h10looper_2pi::getAlpha()","Magnitude of G=%f",TMath::Sqrt(v_G.Dot(v_G)));
	Info("h10looper_2pi::getAlpha()","Angle of G with uv_Gp=%f",ACos(v_G.Dot(uv_Gp))*RadToDeg());*/

	//! step 2.
	float aB=TMath::Sqrt( 1/( 1-TMath::Power(uv_Bf.Dot(uv_Bp),2) ) );
	float bB=-uv_Bf.Dot(uv_Bp)*aB;
	TVector3 v_B=aB*uv_Bf + bB*uv_Bp;
	/*Info("h10looper_2pi::getAlpha()","Magnitude of B=%f",TMath::Sqrt(v_B.Dot(v_B)));
	Info("h10looper_2pi::getAlpha()","Angle of B with uv_Bp=%f",ACos(v_B.Dot(uv_Bp))*RadToDeg());*/

	//! step 3.
	TVector3 v_D=v_G.Cross(v_B);

	//! step 4.
	float alpha=-9999;
	float v_D_angle_uv_Gp=TMath::ACos(v_D.Unit().Dot(uv_Gp))*TMath::RadToDeg();
	//Info("h10looper_2pi::getAlpha()","angle=%f",v_D_angle_uv_Gp);
	(v_D_angle_uv_Gp<90?v_D_angle_uv_Gp=floor(v_D_angle_uv_Gp):v_D_angle_uv_Gp=ceil(v_D_angle_uv_Gp));
	//Info("h10looper_2pi::getAlpha()","angle after floor/ceil=%d",int(v_D_angle_uv_Gp));
	if (int(v_D_angle_uv_Gp)==0){
		alpha=TMath::ACos(v_G.Dot(v_B))*TMath::RadToDeg();
	}else if (int(v_D_angle_uv_Gp)==180){
		alpha=( 2*TMath::Pi()-TMath::ACos(v_G.Dot(v_B)) )*TMath::RadToDeg() ;
	}else{
		Info("h10looper_2pi::getAlpha()","Error: v_D is niether colinear nor anticolinear(angle=%f) with uv_Gp! Setting alpha=-9999",v_D_angle_uv_Gp);
	}

	//! step 5.
	return alpha;
}

void h10looper_2pi::setup_d2pi(){
	TDirectory* dir_d2pi=_fout->mkdir("d2pi");

	//! 1. hq2w
	dir_d2pi->mkdir("Q2W")->cd();
	_hq2w_prec = new TH2F("hq2Vw_prec","Q^{2} vs. W",600,0,3.6,100,0,6); //to match Evgeny;s q2VW bng (Q2/bin, W/bin) = (0.06[GeV^2/bin], 0.006[GeV/bin])
	_hq2w_pstc = new TH2F("hq2Vw_pstc","Q^{2} vs. W",600,0,3.6,100,0,6); //to match Evgeny;s q2VW bng (Q2/bin, W/bin) = (0.06[GeV^2/bin], 0.006[GeV/bin])

	//! 2 Create hmm2ppip, hmmppip,hmm2ppip[NBINS_WCRS],hmmppip[NBINS_WCRS]
	if (_seq=="recon"){
		TDirectory* dir_MM=dir_d2pi->mkdir("MM");
		dir_MM->mkdir("pre_cut")->cd();
		_hmm2_prec_fW=new TH1F("hmm2ppip_prec_fW", "Missing Mass2 of p,#pi^{+}", 100,-0.50,1.00);//0.16));
    	_hmm_prec_fW =new TH1F("hmmppip_prec_fW",  "Missing Mass of p,#pi^{+}",  100,-0.50,1.00);//0.40));
		_hmm2_prec=new TH1F*[NBINS_WCRS];
		_hmm_prec =new TH1F*[NBINS_WCRS];
		for(int iw=0;iw<NBINS_WCRS;iw++){
			_hmm2_prec[iw]=new TH1F(TString::Format("hmm2ppip_prec_%d",iw+1), "Missing Mass2 of p,#pi^{+}", 100,-0.50,1.00);
			_hmm_prec[iw]= new TH1F(TString::Format("hmmppip_prec_%d",iw+1),  "Missing Mass of p,#pi^{+}",  100,-0.50,1.00);
		}
		dir_MM->mkdir("pst_cut")->cd();
		_hmm2_pstc_fW=new TH1F("hmm2ppip_pstc_fW", "Missing Mass2 of p,#pi^{+}", 100,-0.50,1.00);//0.16));
    	_hmm_pstc_fW =new TH1F("hmmppip_pstc_fW",  "Missing Mass of p,#pi^{+}",  100,-0.50,1.00);//0.40));
		_hmm2_pstc=new TH1F*[NBINS_WCRS];
		_hmm_pstc =new TH1F*[NBINS_WCRS];
		for(int iw=0;iw<NBINS_WCRS;iw++){
			_hmm2_pstc[iw]=new TH1F(TString::Format("hmm2ppip_pstc_%d",iw+1), "Missing Mass2 of p,#pi^{+}", 100,-0.50,1.00);
			_hmm_pstc[iw]= new TH1F(TString::Format("hmmppip_pstc_%d",iw+1),  "Missing Mass of p,#pi^{+}",  100,-0.50,1.00);
		}
	}

	//! 3. Create h7[NBINS_WCRS][NVST]
	dir_d2pi->cd();
	Int_t hdim=8;
	struct h7_bng{
    	int bins;
    	float xmin;
    	float xmax;
    } bngQ2, bngW, bngMppip, bngMppim, bngMpippim, bngTheta, bngPhi, bngAlpha;

	//! First set up binning for h7 that is independent of
	//! Wcrs-bng:
	//! +Q2,theta,phi and alpha
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

	//! Now create h7[NBINS_WCRS][NVST] and
	//! in each specific Wcrs-bin, pick appropriate binning for:
	//! + W,M1 and M2
	_h8=new THnSparse**[NBINS_WCRS];
	for(int iw=0;iw<NBINS_WCRS;iw++){
		_h8[iw]=new THnSparse*[NVST];

		//! pick binning for W, M1 and M2
		bngW.bins=WCRSBIN[iw].nbins;
		bngW.xmin=WCRSBIN[iw].xmin;
		bngW.xmax=WCRSBIN[iw].xmax;
			
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
		/* Varset 1 (Delta++) */
		//                    {  h, Q2,         W,         Mppip,         Mpippim,         theta_pim,     phi_pim,     alpha[p'pip][ppim]}
		Int_t bins1[]    =    {  3, bngQ2.bins, bngW.bins, bngMppip.bins, bngMpippim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
		Double_t xmin1[] =    { -1, bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMpippim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
		Double_t xmax1[] =    {  2, bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMpippim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
		_h8[iw][0] = new THnSparseF(TString::Format("h8_%d_%d",iw+1,1), 
		"h, Q^{2}, W, M_{p#pi^{+}}, M_{#pi^{+}#pi^{-}}, #theta_{#pi^{-}}, #phi_{#pi^{-}}, #alpha_{[p^{'}#pi^{+}][p#pi^{-}]}", 
		hdim, bins1, xmin1, xmax1);
		//! If so specified, make variable Q2-binning
		if (_use_Q2_var_binw_bng)_h8[iw][0]->GetAxis(1)->Set(NBINS_Q2_VAR_BINW_BNG,Q2_BINS_VAR_BINW_BNG);
		_h8[iw][0]->Sumw2();
		gDirectory->Append(_h8[iw][0]);
		
		/* Varset 2 (Rho) */
		//                    {  h, Q2,         W,         Mppip,         Mpippim,         theta_p,       phi_p,       alpha[pippim][p'p]}
		Int_t bins2[]    =    {  3, bngQ2.bins, bngW.bins, bngMppip.bins, bngMpippim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
		Double_t xmin2[] =    { -1, bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMpippim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
		Double_t xmax2[] =    {  2, bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMpippim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
		_h8[iw][1] = new THnSparseF(TString::Format("h8_%d_%d",iw+1,2), 
		"h, Q^{2}, W, M_{p#pi^{+}}, M_{#pi^{+}#pi^{-}}, #theta_{p}, #phi_{p}, #alpha_{[#pi^{+}#pi^{-}][pp^{'}]}", 
		hdim, bins2, xmin2, xmax2);
		//! If so specified, make variable Q2-binning
		if (_use_Q2_var_binw_bng)_h8[iw][1]->GetAxis(1)->Set(NBINS_Q2_VAR_BINW_BNG,Q2_BINS_VAR_BINW_BNG);
		_h8[iw][1]->Sumw2();
		gDirectory->Append(_h8[iw][1]);
			
		/* Varset 3 (Delta0) */
		//                    {  h, Q2,         W,         Mppip,         Mppim,         theta_pip,     phi_pip,     alpha[p'pim][ppip]}
		Int_t bins3[]    =    {  3, bngQ2.bins, bngW.bins, bngMppip.bins, bngMppim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
		Double_t xmin3[] =    { -1, bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMppim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
		Double_t xmax3[] =    {  2, bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMppim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
		_h8[iw][2] = new THnSparseF(TString::Format("h8_%d_%d",iw+1,3), 
		"h, Q^{2}, W, M_{p#pi^{+}}, M_{p#pi^{-}}, #theta_{#pi^{+}}, #phi_{#pi^{+}}, #alpha_{[p^{'}#pi^{-}][p#pi^{+}]}", 
		hdim, bins3, xmin3, xmax3);
		//! If so specified, make variable Q2-binning
		if (_use_Q2_var_binw_bng)_h8[iw][2]->GetAxis(1)->Set(NBINS_Q2_VAR_BINW_BNG,Q2_BINS_VAR_BINW_BNG);
		_h8[iw][2]->Sumw2();
		gDirectory->Append(_h8[iw][2]);
	}
	//! Q2,W cut values 
	if (_use_thesis_Q2W){
		_Q2_MIN=1.50,_Q2_MAX=5.00,_W_MIN=1.400,_W_MAX=2.125;
	}else{
		//! [12-07-15] To compare with EI' ana of sim-e16-EI
		_Q2_MIN=2.00,_Q2_MAX=2.40,_W_MIN=1.400,_W_MAX=2.00;
	}
	//!MMcut values
	/*if (_expt=="e1f"){
		_mm2ppip_l=-0.00,_mm2ppip_h=0.04; //0,0.04
	}else if (_expt=="e16"){
		_mm2ppip_l=-0.04,_mm2ppip_h=0.06;
	}*/
	if (_use_MM2_cut_EI){
		_mm2ppip_l=-0.04,_mm2ppip_h=0.06; 
	}else if(_use_MM2_cut_SS){
		_mm2ppip_l=-0.16,_mm2ppip_h=0.16;
	}else{
		_mm2ppip_l=-0.00,_mm2ppip_h=0.04;
	}
}

void h10looper_2pi::fill_h8(){
	int iw=GetCrsWBinIdx(_W);
	
	int h=evthel;
	//! + If _dtyp=="sim", then overwrite h with 0.
	//! + sim-h10 should contain 0, but from what I checked, evthel=7 for sim-h10
	if (_dtyp=="sim"){
		h=0;
	}
	/* Varset 1 (Delta++)*/
	double coord1[] = {h,_Q2,_W,_M_ppip,_M_pippim,_theta_cms_pim,_phi_cms_pim,_alpha_1};
	_h8[iw][0]->Fill(coord1);
	
	/* Varset 2 (Rho)*/
	double coord2[] = {h,_Q2,_W,_M_ppip,_M_pippim,_theta_cms_p,_phi_cms_p,_alpha_2};
	_h8[iw][1]->Fill(coord2);
	
	/* Varset 3 (Delta0)*/
	double coord3[] = {h,_Q2,_W,_M_ppip,_M_ppim,_theta_cms_pip,_phi_cms_pip,_alpha_3};
	_h8[iw][2]->Fill(coord3);

	return;
}

void h10looper_2pi::e16_ST_Loop_with_efid(){
	Info("h10looper_2pi::e16_ST_Loop_with_efid()","");

	if (fChain == 0) return;

	//! output object for EFID
	_fout->mkdir("efid")->cd();
	_hefid=new TH1D("hefid","EFID statistics",NUM_EFID_STATS,0.5,NUM_EFID_STATS+0.5);
	_hefid->GetXaxis()->SetBinLabel(EFID_TOT,"total");
	_hefid->GetXaxis()->SetBinLabel(EFID_IN,"infid");
	//! phi vs. theta hists
	_hefid_e=new TH2F*[2];
	for (int i=0;i<2;i++){
		TString name_sfx;
		if      (i==0) name_sfx="prec";
		else if (i==1) name_sfx="pstc";
		_hefid_e[i]  = new TH2F(TString::Format("h_e_phiVtheta_%s",name_sfx.Data()),  "#phi vs. #theta for e-",  100,0,60,  100,-30,330);
	}

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<_nentries_to_proc;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		if ( (jentry+1)%100000==0 ) {
			Info("h10looper_2pi::e16_ST_Loop_with_efid()","Processing event %llu",jentry+1);
			Info("h10looper_2pi::e16_ST_Loop_with_efid()","Fraction of events processed %.2f%%",((float)(jentry+1)/_nentries_to_proc)*100);
		}

		//! Reconcile h10 Branch binding for e16:exp
		//! or
		//! if _do_reconcile specified by user (for example for sim-e16-EI)
		if (_do_reconcile || (_expt=="e16" && _dtyp=="exp")) {
			Reconcile();
		}

		//! Begin processing event
		reset_ekin();
		reset_hkin();
		_hevt->Fill(EVT_TRG);
		set_ekin();
		set_hkin();
			
		//! Event selection
		//! EFID
		_hefid->Fill(EFID_TOT);
		if (e16_ST_pass_efid()){
			_hefid->Fill(EFID_IN);
			_hevt->Fill(EVT_E_INFID);
		
			//! Q2-W kinematic cut
			_hq2w_prec->Fill(_W,_Q2);
			if (!(_Q2>=_Q2_MIN && _Q2<_Q2_MAX && _W>_W_MIN && _W<_W_MAX)) continue;
			_hevt->Fill(EVT_Q2W_KIN_PASS);
			_hq2w_pstc->Fill(_W,_Q2);

			_hevt->Fill(EVT_2PI);
			fill_h8();
		}
	}//end event loop

	//! Write _fout
	_fout->Write();
}//end Loop()

/* Currently implemented for top2'*/
bool h10looper_2pi::pass_eff(){
	//Info("h10looper_2pi::pass_eff()","");
	bool ret=kFALSE;
	bool e_ineff=kFALSE;
	bool p_ineff=kFALSE;
	bool pip_ineff=kFALSE;

	//Info("h10looper_2pi::pass_eff()","sectors:%d, %d, %d",_sector_e,_sector_p,_sector_pip);

	//!prec
	_hthetavp[0][_sector_e-1][0]->Fill(_p_e,_theta_e);
	_hthetavp[1][_sector_p-1][0]->Fill(_p_p,_theta_p);
	_hthetavp[2][_sector_pip-1][0]->Fill(_p_pip,_theta_pip);

	//! check if e- ineff
	if (pass_theta_vs_p("e")){
		_heff->Fill(EFF_E_PASS);
		e_ineff=kTRUE;
		_hthetavp[0][_sector_e-1][1]->Fill(_p_e,_theta_e);
	}
	//! check if p ineff
	if (pass_theta_vs_p("p")){
		_heff->Fill(EFF_P_PASS);
		p_ineff=kTRUE;
		_hthetavp[1][_sector_p-1][1]->Fill(_p_p,_theta_p);
	}
	//! check if pip ineff
	if (pass_theta_vs_p("pip")){
		_heff->Fill(EFF_PIP_PASS);
		pip_ineff=kTRUE;
		_hthetavp[2][_sector_pip-1][1]->Fill(_p_pip,_theta_pip);
	}

	ret=e_ineff && p_ineff && pip_ineff;
	return ret;
}

/* Currently implemented for top2'*/
bool h10looper_2pi::pass_scpd(){
	//Info("h10looper_2pi::pass_scpd()","");
	bool ret=kFALSE;
	bool e_inscpd=kFALSE;
	bool p_inscpd=kFALSE;
	bool pip_inscpd=kFALSE;

	//Info("h10looper_2pi::pass_scpd()","sectors:%d, %d, %d",_sector_e,_sector_p,_sector_pip);

	//!prec
	_hpdl[0][0]->Fill(_sc_pd_e);
	_hpdl[1][0]->Fill(_sc_pd_p);
	_hpdl[2][0]->Fill(_sc_pd_pip);
	_hthetavp2[0][_sector_e-1][0]->Fill(_p_e,_theta_e);
	_hthetavp2[1][_sector_p-1][0]->Fill(_p_p,_theta_p);
	_hthetavp2[2][_sector_pip-1][0]->Fill(_p_pip,_theta_pip);

	//! check if e- inscpd
	//Info("h10looper_2pi::pass_scpd()","sc_pd_e=%d",_sc_pd_e);
	if (!is_scpd_bad("e")){
		_hscpd->Fill(SCPD_E_PASS);
		e_inscpd=kTRUE;
		_hpdl[0][1]->Fill(_sc_pd_e);
		_hthetavp2[0][_sector_e-1][1]->Fill(_p_e,_theta_e);
	}
	//! check if p inscpd
	//Info("h10looper_2pi::pass_scpd()","sc_pd_p=%d",_sc_pd_p);
	if (!is_scpd_bad("p")){
		_hscpd->Fill(SCPD_P_PASS);
		p_inscpd=kTRUE;
		_hpdl[1][1]->Fill(_sc_pd_p);
		_hthetavp2[1][_sector_p-1][1]->Fill(_p_p,_theta_p);
	}
	//! check if pip inscpd
	//Info("h10looper_2pi::pass_scpd()","sc_pd_pip=%d",_sc_pd_pip);
	if (!is_scpd_bad("pip")){
		_hscpd->Fill(SCPD_PIP_PASS);
		pip_inscpd=kTRUE;
		_hpdl[2][1]->Fill(_sc_pd_pip);
		_hthetavp2[2][_sector_pip-1][1]->Fill(_p_pip,_theta_pip);
	}

	ret=e_inscpd && p_inscpd && pip_inscpd;
	return ret;
}

/*
[02-28-16] hack-pcorr for '_make_h10_skim_SS'
*/
void h10looper_2pi::set_hkin_use_passed_lv(TLorentzVector lvP, TLorentzVector lvPip,int h10idx_p/*=-1*/, int h10idx_pip/*=-1*/){
	_lvP=lvP;
	_lvPip=lvPip;
	//! Now set _lvPim
	_lvPim = (_lvW-(_lvP+_lvPip));

	//! Set up Lab frame components of kin.
	_p_p=_lvP.P();
	_p_pip=_lvPip.P();
	_p_pim=_lvPim.P();
	_theta_p=_lvP.Theta()*TMath::RadToDeg();
	_theta_pip=_lvPip.Theta()*TMath::RadToDeg();
	_theta_pim=_lvPim.Theta()*TMath::RadToDeg();
	//! phi: 1st get in range [-180,180]
	float phi_p=_lvP.Phi()*TMath::RadToDeg();// [-180,180]
	float phi_pip=_lvPip.Phi()*TMath::RadToDeg();// [-180,180]
	float phi_pim=_lvPim.Phi()*TMath::RadToDeg();// [-180,180]
	//! phi: now convert to [-30,330]
	_phi_p=phi_p<-30?phi_p+360:phi_p; // [-30,330]
	_phi_pip=phi_pip<-30?phi_pip+360:phi_pip; // [-30,330]
	_phi_pim=phi_pim<-30?phi_pim+360:phi_pim; // [-30,330]
	//! sector
	_sector_p=get_sector(_phi_p);
	_sector_pip=get_sector(_phi_pip);
	_sector_pim=get_sector(_phi_pim);
	//! paddle
	if (_seq=="recon"){
		_sc_pd_p=100*dc_sect[dc[h10idx_p]-1]+sc_pd[sc[h10idx_p]-1];
		_sc_pd_pip=100*dc_sect[dc[h10idx_pip]-1]+sc_pd[sc[h10idx_pip]-1];
		//!_sc_pd_pim= thus far, not identifying pim
	}

	//! + The following code taken directly from ProcD2pi::UpdateD2pi
	//!   and therein the part to calculate Varset kinematics
	//! + The only changes relate to certain variable names that have 
	//!   changed and code readability etc.

	//! Varsets kinematics
	//! Calculate rotation: taken from Evan's phys-ana-omega on 08-05-13
	TVector3 uz = _lvQ.Vect().Unit();
    TVector3 ux = (_lvE0.Vect().Cross(_lvE1.Vect())).Unit();
    ux.Rotate(-TMath::Pi()/2,uz);
    TRotation r3;// = new TRotation();
    r3.SetZAxis(uz,ux).Invert();
    //! _w and _q are in z-direction
    TVector3 boost(-1*_lvW.BoostVector());
    TLorentzRotation r4(r3); //*_boost);
    r4 *= boost; //*_3rot;
	
	_lvQCMS   = _lvQ;
	_lvP0CMS  = _lvP0;
	_lvPCMS   = _lvP;
	_lvPipCMS = _lvPip;
	_lvPimCMS = _lvPim;
	//! Rotate and Boost
	_lvQCMS.Transform(r4);
	_lvP0CMS.Transform(r4);
	_lvPCMS.Transform(r4);
	_lvPipCMS.Transform(r4);
	_lvPimCMS.Transform(r4);
	
	_M_ppip = (_lvPCMS + _lvPipCMS).Mag();
	_M_ppim = (_lvPCMS + _lvPimCMS).Mag();
	_M_pippim = (_lvPipCMS + _lvPimCMS).Mag();
	
	//! Should be directly be able to get the Theta() angle since
	//! virtual photon defines the coordinate system
	_theta_cms_p=_lvPCMS.Theta()*TMath::RadToDeg();
	_theta_cms_pip=_lvPipCMS.Theta()*TMath::RadToDeg();
	_theta_cms_pim=_lvPimCMS.Theta()*TMath::RadToDeg();
	
	_phi_cms_p=getPhi(_lvPCMS);
	_phi_cms_pip=getPhi(_lvPipCMS);
	_phi_cms_pim=getPhi(_lvPimCMS);

	//! alpha angle
	//! Following vectors are used in the calculation of alpha angle; see doc. for getAlpha()
	TVector3 G_f(0,0,0);
	TVector3 G_p(0,0,0);
	TVector3 B_f(0,0,0);
	TVector3 B_p(0,0,0);


	G_f=uz*-1; //! Gamma vector used in calculation of alpha always "follows" -uz

	//! alpha[p',pip][p,pim]
	G_p=_lvPimCMS.Vect().Unit();
	B_f=_lvPipCMS.Vect().Unit();
	B_p=G_p;
	_alpha_1=getAlpha(G_f,G_p,B_f,B_p);
	//! alpha[pip,pim][p,p']
	G_p=_lvPCMS.Vect().Unit();
	B_f=_lvPipCMS.Vect().Unit();
	B_p=G_p;
	_alpha_2=getAlpha(G_f,G_p,B_f,B_p);
	//! alpha[p',pim][p,pip]
	G_p=_lvPipCMS.Vect().Unit();
	B_f=_lvPimCMS.Vect().Unit();//_lvPCMS.Vect().Unit();
	B_p=G_p;
	_alpha_3=getAlpha(G_f,G_p,B_f,B_p);

	//! h10idx
	_h10idx_p=h10idx_p;
	_h10idx_pip=h10idx_pip;
	//_h10idx_pim=TBD;
}
/*
[02-28-16] End: hack-pcorr for '_make_h10_skim_SS'
*/