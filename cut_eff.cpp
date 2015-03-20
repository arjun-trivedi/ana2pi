#include "cut_eff.h"
#include "particle_constants.h"
#include <TCut.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TTree.h>
/*#include <TF1.h>
#include <TString.h>*/

#include <cstdio>
/*#include <iostream>
using namespace std;*/

using namespace ParticleConstants;

CutEff::CutEff():_prtcl{"e","p","pip","pim"},
_theta_min{0,0,0,0},_theta_max{60,60,120,120},
_p_min{1,0,0,0},_p_max{5,4,3,3}
{
  setup_cuts();
}

CutEff::~CutEff()
{
  delete _cut_lw,_cut_hg;
}

Bool_t CutEff::InEfficientRegion(Int_t pid,Int_t sector,Float_t theta,Float_t p)
{
  //! + Assume that particle will be in efficient region
  //! + The ret value will be set to false if the functions finds out that
  //!   the particle is in the inefficient region of the detector
  Bool_t ret=kTRUE; 
  //! prtcl index
  Int_t i=-1; 
  if (pid==ELECTRON) i=0;
  if (pid==PROTON) i=1;
  if (pid==PIP) i=2;
  if (pid==PIM) i=3;
  //! sctr index
  Int_t isctr=sector-1;

  //! See if prtcl is in detector's inefficient region
  //! + e,p and pim have only 1 inefficient region per sector
  //! + pip has more than 1 inefficient region(sectors 2,3 and 6)
  if (_prtcl[i]=="e"||_prtcl[i]=="p"||_prtcl[i]=="pim"){
    if (_cut_lw[i][isctr]!=NULL && _cut_hg[i][isctr]!=NULL){
      Float_t theta_lw=_cut_lw[i][isctr][0]->Eval(p);
      Float_t theta_hg=_cut_hg[i][isctr][0]->Eval(p);
      if ( (theta>=theta_lw)&&(theta<=theta_hg) ){//! prtcl is in inefficient region
        ret=kFALSE;
      }
    }
  }else if (_prtcl[i]=="pip"){
    if (_cut_lw[i][isctr]!=NULL && _cut_hg[i][isctr]!=NULL){
      Int_t ncuts=1;
      if (isctr+1==2) ncuts=2;
      if (isctr+1==3) ncuts=4;
      if (isctr+1==6) ncuts=2;
      for (int j=0;j<ncuts;j++){
          Float_t theta_lw=_cut_lw[i][isctr][j]->Eval(p);
          Float_t theta_hg=_cut_hg[i][isctr][j]->Eval(p);
          if ( (theta>=theta_lw)&&(theta<=theta_hg) ){//! prtcl is in inefficient region
            ret=kFALSE;
            break; //break out of the for loop is in any of sector's inefficient regions
          }
      }
    }
  }
  return ret;
}

void CutEff::setup_cuts(){
  for (int i=0;i<4;i++){
    Double_t pmin=_p_min[i];
    Double_t pmax=_p_max[i];
    Int_t thetamin=_theta_min[i];
    Int_t thetamax=_theta_max[i];
    //printf("pmin,pax,thetamin,thetamax=%.1f,%.1f,%d,%d\n",pmin,pmax,thetamin,thetamax);
    for (int isctr=0;isctr<6;isctr++){
      TString name=TString::Format("%s_%d",_prtcl[i].Data(),isctr+1).Data();
      if (_prtcl[i]=="e"){
        if (isctr+1==1||isctr+1==2||isctr+1==5||isctr+1==6){//! no cuts
          _cut_lw[i][isctr]=NULL;
          _cut_hg[i][isctr]=NULL;		
        }else if (isctr+1==3){
          _cut_lw[i][isctr]=new TF1*[1];
          _cut_hg[i][isctr]=new TF1*[1];          
          _cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"30-1*x",2.0,3.5);
          _cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"32-1*x",2.0,3.5);        
        }else if (isctr+1==4){
          _cut_lw[i][isctr]=new TF1*[1];
          _cut_hg[i][isctr]=new TF1*[1];
          _cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"21.0-1*x",2.0,4.5);
          _cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"24.5-1*x",2.0,4.5);
        }
      }
      if (_prtcl[i]=="p"){
        if (isctr+1==1||isctr+1==4||isctr+1==6){//! no cuts
          _cut_lw[i][isctr]=NULL;
          _cut_hg[i][isctr]=NULL;
        }else if (isctr+1==2){
          _cut_lw[i][isctr]=new TF1*[1];
          _cut_hg[i][isctr]=new TF1*[1];
          _cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"16+8*x-1*x*x",0.5,3.0);
          _cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"19+8*x-1*x*x",0.5,3.0); 
        }else if (isctr+1==3){
          _cut_lw[i][isctr]=new TF1*[1];
          _cut_hg[i][isctr]=new TF1*[1];
          _cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"5+8*x-1*x*x",0.5,3.5);
          _cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"8+8*x-1*x*x",0.5,3.5);
        }else if (isctr+1==5){
          _cut_lw[i][isctr]=new TF1*[1];
          _cut_hg[i][isctr]=new TF1*[1];
          _cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"19+5*x-1*x*x",0.5,3.0);
          _cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"22+5*x-1*x*x",0.5,3.0);
        }
      }
      if (_prtcl[i]=="pip"){
        if (isctr+1==1||isctr+1==4||isctr+1==5){//! no cuts
          _cut_lw[i][isctr]=NULL;
          _cut_hg[i][isctr]=NULL;
        }else if (isctr+1==2){
          _cut_lw[i][isctr]=new TF1*[2];
          _cut_hg[i][isctr]=new TF1*[2];
          _cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"50+105*x-90*x*x",0.03,0.7);
          _cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),TString::Format("%d",thetamax),0.03,0.7);
          _cut_lw[i][isctr][1]=new TF1(TString::Format("%s_%d_lw",name.Data(),2),"-13+37*x-5*x*x",0.5,2.5);
          _cut_hg[i][isctr][1]=new TF1(TString::Format("%s_%d_hg",name.Data(),2),"-7+38*x-4*x*x",0.4,2.5);
        }else if (isctr+1==3){
          _cut_lw[i][isctr]=new TF1*[4];
          _cut_hg[i][isctr]=new TF1*[4];
          _cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"80+37*x-20*x*x",0.0,0.7);
          _cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"85+40*x-15*x*x",0.0,0.7);
          _cut_lw[i][isctr][1]=new TF1(TString::Format("%s_%d_lw",name.Data(),2),"64.5+32*x-20*x*x",0.0,0.7);
          _cut_hg[i][isctr][1]=new TF1(TString::Format("%s_%d_hg",name.Data(),2),"68+35*x-18*x*x",0.0,0.7);
          _cut_lw[i][isctr][2]=new TF1(TString::Format("%s_%d_lw",name.Data(),3),"15+30*x-7*x*x",0.1,1.5);
          _cut_hg[i][isctr][2]=new TF1(TString::Format("%s_%d_hg",name.Data(),3),"20.5+30*x-7*x*x",0.1,1.5);
          _cut_lw[i][isctr][3]=new TF1(TString::Format("%s_%d_lw",name.Data(),4),"-13+25*x-5*x*x",0.6,2.5);
          _cut_hg[i][isctr][3]=new TF1(TString::Format("%s_%d_hg",name.Data(),4),"-8.0+25*x-5*x*x",0.6,2.5);
        }else if (isctr+1==6){
          _cut_lw[i][isctr]=new TF1*[2];
          _cut_hg[i][isctr]=new TF1*[2];
          _cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"64.5+32*x-20*x*x",0.0,0.7);
          _cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"68+35*x-18*x*x",0.0,0.7);
          _cut_lw[i][isctr][1]=new TF1(TString::Format("%s_%d_lw",name.Data(),2),"0+30*x-7*x*x",0.1,2.0);
          _cut_hg[i][isctr][1]=new TF1(TString::Format("%s_%d_hg",name.Data(),2),"12+30*x-7*x*x",0.1,2.0);
        }
      }
      if (_prtcl[i]=="pim"){
        if (isctr+1==1||isctr+1==4||isctr+1==5||isctr+1==6){//! no cuts
          _cut_lw[i][isctr]=NULL;
          _cut_hg[i][isctr]=NULL;
        }else if (isctr+1==2){
          _cut_lw[i][isctr]=new TF1*[1];
          _cut_hg[i][isctr]=new TF1*[1];
          _cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"55-25*x+5*x*x",0.2,2);
          _cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"62-25*x+5*x*x",0.2,2);
        }else if (isctr+1==3){
          _cut_lw[i][isctr]=new TF1*[1];
          _cut_hg[i][isctr]=new TF1*[1];
          _cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"49-25*x+8*x*x",0.2,2);
          _cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"51.5-25*x+8*x*x",0.2,2);
        }
      }
    }
  }
}