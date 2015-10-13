#define h10looper_2pi_cxx
#include "h10looper_2pi.h"

#include "wrpr_cut_fid_e16.h"
#include "h7_bng.h"

#include <TLorentzRotation.h>

h10looper_2pi::h10looper_2pi(TString h10type,TChain* h10chain,TString fout_name,Long64_t nentries,
						     TString adtnl_cut_opt) : h10looper_e1f(h10type,h10chain,fout_name,nentries,adtnl_cut_opt)
{
	Info("h10looper_2pi::h10looper_2pi","Setting up h10looper_2pi...\n");

	//! output objects
	if (_seq=="recon"){
		_hpfid=new TH1D("hpfid","PFID statistics",NUM_PFID_STATS,0.5,NUM_PFID_STATS+0.5);
		_hpfid->GetXaxis()->SetBinLabel(PFID_TOT,"total");
		_hpfid->GetXaxis()->SetBinLabel(PFID_P_IN,"p infid");
		_hpfid->GetXaxis()->SetBinLabel(PFID_PIP_IN,"pip infid");
		_hpfid->GetXaxis()->SetBinLabel(PFID_P_AND_PIP_IN,"p+pip infid");
	}

	//! d2pi
	setup_d2pi();

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
	Info("","---------------------------------------");
	Info("h10looper_2pi::h10looper_2pi","Done setting up h10looper_2pi\n");
	
}

h10looper_2pi::~h10looper_2pi()
{
	Info("h10looper_2pi::~h10looper_2pi","");
	
	delete _hpfid;

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
	delete[] _h7;

	Info("h10looper_2pi::~h10looper_2pi","Done.");
}

void h10looper_2pi::Loop(){
	Info("h10looper_2pi::Loop()","");

	if (fChain == 0) return;

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<_nentries_to_proc;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		if ( (jentry+1)%100000==0 ) {
			Info("h10looper_2pi::Loop()","Processing event %llu",jentry+1);
			Info("h10looper_2pi::Loop()","Fraction of events processed %.2f%%",((float)(jentry+1)/_nentries_to_proc)*100);
		}
		//! Begin processing event
		reset_ekin();
		reset_hkin();
		if (_seq=="recon"){
			_hevt->Fill(EVT_TRG);
			//! EID
			if (evt_trigger_electron()){
				_hevt->Fill(EVT_E);
				set_ekin();
				//! EFID
				if (electron_infid()){
					_hevt->Fill(EVT_E_INFID);
					//! pcorr
					if (_dtyp=="exp"){
						mom_corr_electron();
						set_ekin();
					}
					//! PID
					//! only top2':
					//!   + only p and pip need to be detected
					//!   + pim is always reconstructed using lvPim=lvW-(lvP+lvPip)
					//!     i.e. it does not matter if it is detected or missing
					int h10idx_p=found_hadron("p");
					int h10idx_pip=found_hadron("pip");
					if (h10idx_p>=1 && h10idx_pip>=1){ 
						_hpid->Fill(PID_P_AND_PIP_FOUND);
						_hevt->Fill(EVT_P_PIP);
						set_hkin(h10idx_p, h10idx_pip);
						//! PFID
						if (proton_infid() && pip_infid()){
							_hpfid->Fill(PFID_P_AND_PIP_IN);
							_hevt->Fill(EVT_P_PIP_INFID);
							
							//! Event selection

							//! Q2-W kinematic cut
							_hq2w_prec->Fill(_W,_Q2);
							if (!(_Q2>=1.25 && _Q2<5.25 && _W>1.300 && _W<2.125)) continue;
							_hevt->Fill(EVT_Q2W_KIN_PASS);
							_hq2w_pstc->Fill(_W,_Q2);

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

							//! Apply MM cut
							if (mm2ppip>0 && mm2ppip<0.04){
								_hevt->Fill(EVT_2PI);
								_hmm2_pstc_fW->Fill(mm2ppip);
								_hmm_pstc_fW->Fill(mmppip);
								_hmm2_prec[iw]->Fill(mm2ppip);
								_hmm_prec[iw]->Fill(mmppip);
								fill_h7();
							}
						}
					}
					
				}
			}
		}else if (_seq=="thrown"){
			_hevt->Fill(EVT_TRG);
			set_ekin();
			set_hkin();
			
			//! Event selection

			//! Q2-W kinematic cut
			_hq2w_prec->Fill(_W,_Q2);
			if (!(_Q2>=1.25 && _Q2<5.25 && _W>1.300 && _W<2.125)) continue;
			_hevt->Fill(EVT_Q2W_KIN_PASS);
			_hq2w_pstc->Fill(_W,_Q2);

			_hevt->Fill(EVT_2PI);
			fill_h7();
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
    //! Varset kinematics
    _M_ppip=_M_ppim=_M_pippim=0;
    _p_cms_p=  _theta_cms_p=  _phi_cms_p=0;
    _p_cms_pip=_theta_cms_pip=_phi_cms_pip=0;
    _p_cms_pim=_theta_cms_pim=_phi_cms_pim=0;
    _alpha_1=_alpha_2=_alpha_3=0;
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

}

bool h10looper_2pi::proton_infid(){
	_hpfid->Fill(PFID_TOT);
	int sctr_p=get_sector(_phi_p);
	bool ret=Fiducial_e16_hdrn(_theta_p,_phi_p,sctr_p);
	if (ret==kTRUE){
		_hpfid->Fill(PFID_P_IN);	
	} 
	return ret;
}

bool h10looper_2pi::pip_infid(){
	_hpfid->Fill(PFID_TOT);
	int sctr_pip=get_sector(_phi_pip);
	bool ret=Fiducial_e16_hdrn(_theta_pip,_phi_pip,sctr_pip);
	if (ret==kTRUE){
		 _hpfid->Fill(PFID_PIP_IN);
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
		_hmm2_prec_fW=new TH1F("hmm2ppip_prec_fW", "Missing Mass2 of p,#pi^{+}", 100,-0.02,1.00);//0.16));
    	_hmm_prec_fW =new TH1F("hmmppip_prec_fW",  "Missing Mass of p,#pi^{+}",  100, 0.00,1.00);//0.40));
		_hmm2_prec=new TH1F*[NBINS_WCRS];
		_hmm_prec =new TH1F*[NBINS_WCRS];
		for(int iw=0;iw<NBINS_WCRS;iw++){
			_hmm2_prec[iw]=new TH1F(TString::Format("hmm2ppip_prec_%d",iw+1), "Missing Mass2 of p,#pi^{+}", 100,-0.02,1.00);
			_hmm_prec[iw]= new TH1F(TString::Format("hmmppip_prec_%d",iw+1),  "Missing Mass of p,#pi^{+}",  100, 0.00,1.00);
		}
		dir_MM->mkdir("pst_cut")->cd();
		_hmm2_pstc_fW=new TH1F("hmm2ppip_pstc_fW", "Missing Mass2 of p,#pi^{+}", 100,-0.02,1.00);//0.16));
    	_hmm_pstc_fW =new TH1F("hmmppip_pstc_fW",  "Missing Mass of p,#pi^{+}",  100, 0.00,1.00);//0.40));
		_hmm2_pstc=new TH1F*[NBINS_WCRS];
		_hmm_pstc =new TH1F*[NBINS_WCRS];
		for(int iw=0;iw<NBINS_WCRS;iw++){
			_hmm2_pstc[iw]=new TH1F(TString::Format("hmm2ppip_pstc_%d",iw+1), "Missing Mass2 of p,#pi^{+}", 100,-0.02,1.00);
			_hmm_pstc[iw]= new TH1F(TString::Format("hmmppip_pstc_%d",iw+1),  "Missing Mass of p,#pi^{+}",  100, 0.00,1.00);
		}
	}

	//! 3. Create h7[NBINS_WCRS][NVST]
	dir_d2pi->cd();
	Int_t hdim=7;
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
	_h7=new THnSparse**[NBINS_WCRS];
	for(int iw=0;iw<NBINS_WCRS;iw++){
		_h7[iw]=new THnSparse*[NVST];

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
		//                    {  Q2,         W,         Mppip,         Mpippim,         theta_pim,     phi_pim,     alpha[p'pip][ppim]}
		Int_t bins1[]    =    {  bngQ2.bins, bngW.bins, bngMppip.bins, bngMpippim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
		Double_t xmin1[] =    {  bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMpippim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
		Double_t xmax1[] =    {  bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMpippim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
		_h7[iw][0] = new THnSparseF(TString::Format("h7_%d_%d",iw+1,1), 
		"Q^{2}, W, M_{p#pi^{+}}, M_{#pi^{+}#pi^{-}}, #theta_{#pi^{-}}, #phi_{#pi^{-}}, #alpha_{[p^{'}#pi^{+}][p#pi^{-}]}", 
		hdim, bins1, xmin1, xmax1);
		_h7[iw][0]->Sumw2();
		gDirectory->Append(_h7[iw][0]);
		
		/* Varset 2 (Rho) */
		//                    {  Q2,         W,         Mppip,         Mpippim,         theta_p,       phi_p,       alpha[pippim][p'p]}
		Int_t bins2[]    =    {  bngQ2.bins, bngW.bins, bngMppip.bins, bngMpippim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
		Double_t xmin2[] =    {  bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMpippim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
		Double_t xmax2[] =    {  bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMpippim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
		_h7[iw][1] = new THnSparseF(TString::Format("h7_%d_%d",iw+1,2), 
		"Q^{2}, W, M_{p#pi^{+}}, M_{#pi^{+}#pi^{-}}, #theta_{p}, #phi_{p}, #alpha_{[#pi^{+}#pi^{-}][pp^{'}]}", 
		hdim, bins2, xmin2, xmax2);
		_h7[iw][1]->Sumw2();
		gDirectory->Append(_h7[iw][1]);
			
		/* Varset 3 (Delta0) */
		//                    {  Q2,         W,         Mppip,         Mppim,         theta_pip,     phi_pip,     alpha[p'pim][ppip]}
		Int_t bins3[]    =    {  bngQ2.bins, bngW.bins, bngMppip.bins, bngMppim.bins, bngTheta.bins, bngPhi.bins, bngAlpha.bins };
		Double_t xmin3[] =    {  bngQ2.xmin, bngW.xmin, bngMppip.xmin, bngMppim.xmin, bngTheta.xmin, bngPhi.xmin, bngAlpha.xmin };
		Double_t xmax3[] =    {  bngQ2.xmax, bngW.xmax, bngMppip.xmax, bngMppim.xmax, bngTheta.xmax, bngPhi.xmax, bngAlpha.xmax };
		_h7[iw][2] = new THnSparseF(TString::Format("h7_%d_%d",iw+1,3), 
		"Q^{2}, W, M_{p#pi^{+}}, M_{p#pi^{-}}, #theta_{#pi^{+}}, #phi_{#pi^{+}}, #alpha_{[p^{'}#pi^{-}][p#pi^{+}]}", 
		hdim, bins3, xmin3, xmax3);
		_h7[iw][2]->Sumw2();
		gDirectory->Append(_h7[iw][2]);
	}
}

void h10looper_2pi::fill_h7(){
	int iw=GetCrsWBinIdx(_W);
	
	/* Varset 1 (Delta++)*/
	double coord1[] = {_Q2,_W,_M_ppip,_M_pippim,_theta_cms_pim,_phi_cms_pim,_alpha_1};
	_h7[iw][0]->Fill(coord1);
	
	/* Varset 2 (Rho)*/
	double coord2[] = {_Q2,_W,_M_ppip,_M_pippim,_theta_cms_p,_phi_cms_p,_alpha_2};
	_h7[iw][1]->Fill(coord2);
	
	/* Varset 3 (Delta0)*/
	double coord3[] = {_Q2,_W,_M_ppip,_M_ppim,_theta_cms_pip,_phi_cms_pip,_alpha_3};
	_h7[iw][2]->Fill(coord3);

	return;
}