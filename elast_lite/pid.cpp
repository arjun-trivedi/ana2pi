#include "pid.h"
#include <TString.h>
//#include "ep_processor.h"
#include <TChain.h>
#include <TFileCollection.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRegexp.h>
//! Following changes made to pid.cpp in elast_lite version
/*#include "data_h10.h"
#include "data_ana.h"
#include "h10looper.h"*/
#include <TLorentzVector.h> //! Makes Info() work

Pid::Pid(TString dtyp)
{
	Info("Pid::Pid()", "");
	_dtyp=dtyp;
	Info("Pid::Pid()", "dtyp=%s",_dtyp.Data());

	//! Initialize data members
	_prtcl_name=new TString[3];
	_prtcl_name[0]="p";
	_prtcl_name[1]="pip";
    _prtcl_name[2]="pim";
    for (int i=0;i<_NPRTCLS;i++){
    	Info("Pid::Pid()", "particle[%d]=%s",i,_prtcl_name[i].Data());
    }

    _fcut_l=new TF1* [_NPRTCLS];
    _fcut_h=new TF1* [_NPRTCLS];
    for (int i=0;i<_NPRTCLS;i++){
    	TString name;

    	name=TString::Format("fcut_l_%s",_prtcl_name[i].Data());
    	_fcut_l[i]=new TF1(name.Data(),"pol3",0,5);

    	name=TString::Format("fcut_h_%s",_prtcl_name[i].Data());
    	_fcut_h[i]=new TF1(name.Data(),"pol3",0,5);
    }
    
	//! Get cut parameters	
	TString work_space=getenv("WORKSPACE");
	TString fpath=TString::Format("%s/ana2pi/sub_studies/study_pid/dtVp_cuts/%s",work_space.Data(),_dtyp.Data());
	Info("Pid::Pid()", "cut file path=%s",fpath.Data());

	for (int i=0;i<3;i++){
		Info("Pid::Pid()", "particle[%d]=%s",i,_prtcl_name[i].Data());
		TString fname=TString::Format("%s/%s/cutpars.txt",fpath.Data(),_prtcl_name[i].Data());
		Info("Pid::Pid()", "Reading file=%s",fname.Data());
		FILE *fp = fopen(fname.Data(), "r");
		if (fp) {
			char str[200];
			char dscrptn[20];
			float p0,p1,p2,p3;
			while (fgets(str,sizeof(str),fp)) {
				strtok(str, "\n"); //strip off '\n' from str
				Info("Pid::Pid()", "cut pars directly read from file=%s",str);
				sscanf(str, "%s %f %f %f %f", dscrptn,&p0, &p1,&p2,&p3);
				//Info("Pid::Pid()", "dscrptn,&p0,&p1,&p2,&p3 = %s %f %f %f %f",dscrptn,p0,p1,p2,p3);
				if (strstr(dscrptn,"cut_l")!=NULL){
					Info("Pid::Pid()","Initializing cut_l for %s",_prtcl_name[i].Data());
					_fcut_l[i]->SetParameter(0,p0);
					_fcut_l[i]->SetParameter(1,p1);
					_fcut_l[i]->SetParameter(2,p2);
					_fcut_l[i]->SetParameter(3,p3);
					Info("Pid::Pid()", "cut pars= %f %f %f %f",_fcut_l[i]->GetParameter(0),_fcut_l[i]->GetParameter(1),
						                                       _fcut_l[i]->GetParameter(2),_fcut_l[i]->GetParameter(3));
				}else if (strstr(dscrptn,"cut_h")!=NULL){
					Info("Pid::Pid()","Initializing cut_h for %s",_prtcl_name[i].Data());
					_fcut_h[i]->SetParameter(0,p0);
					_fcut_h[i]->SetParameter(1,p1);
					_fcut_h[i]->SetParameter(2,p2);
					_fcut_h[i]->SetParameter(3,p3);
					Info("Pid::Pid()", "cut pars= %f %f %f %f",_fcut_h[i]->GetParameter(0),_fcut_h[i]->GetParameter(1),
						                                       _fcut_h[i]->GetParameter(2),_fcut_h[i]->GetParameter(3));
				}else{
					Info("Pid::Pid()","dscrptn not correct!");
				}
			}
 		}else{
			Info("Pid::Pid()","File cannot be read!");
		}
		fclose(fp);
	}
}

Pid::~Pid()
{
	delete[] _prtcl_name;
	delete[] _fcut_l;
    delete[] _fcut_h;
}

bool Pid::is_proton(float dt,float p){
	float min=_fcut_l[0]->Eval(p);
	float max=_fcut_h[0]->Eval(p);
	/*cout<<"dt="<<dt<<endl;
	cout <<"min,max="<<min<<","max<<endl;*/
	return (dt>=min && dt<=max);
}

bool Pid::is_pip(float dt,float p){
	float min=_fcut_l[1]->Eval(p);
	float max=_fcut_h[1]->Eval(p);
	return (dt>=min && dt<=max);
}

bool Pid::is_pim(float dt,float p){
	float min=_fcut_l[2]->Eval(p);
	float max=_fcut_h[2]->Eval(p);
	return (dt>=min && dt<=max);
}

void Pid::plot_dtVp_cuts(){
	TCanvas* c=new TCanvas(TString::Format("c_%s",_dtyp.Data()),TString::Format("c_%s",_dtyp.Data()));
	c->Divide(2,2);

	for (int i=0;i<_NPRTCLS;i++){
    	TF1* fcutl=_fcut_l[i];
    	TF1* fcuth=_fcut_h[i];
    	fcutl->SetLineColor(kBlue);
    	fcuth->SetLineColor(kRed);
    	fcutl->SetMinimum(-5);
    	fcutl->SetMaximum(5);
    	fcuth->SetMinimum(-5);
    	fcuth->SetMaximum(5);
    	c->cd(i+1);
    	fcutl->Draw();
    	fcuth->Draw("same");
    	TLegend* l=new TLegend(0.1,0.7,0.48,0.9);
		l->AddEntry(fcutl->GetName(),TString::Format("%s_%s",fcutl->GetName(),_dtyp.Data()),"l");
		l->AddEntry(fcuth->GetName(),TString::Format("%s_%s",fcuth->GetName(),_dtyp.Data()),"l");
		l->Draw();
    	c->SaveAs(TString::Format("/tmp/%s.eps",_dtyp.Data()));
    }
    //delete c; 
}
