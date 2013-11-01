#include "anahND.h"
#include <TFile.h>
#include <TKey.h>

#include <TRandom.h>
#include <TStyle.h>

#include <cstdio>
#include <fstream>


//a: acceptance study 
//d: dead-cells study
Int_t rNum;
void anahND(char* fname, TString type="a", Bool_t mDraw=kFALSE){
	_type = type;
	gStyle->SetOptStat("nMReiou");
	gStyle->SetErrorX(0);
 
	TFile *f_exp = NULL;
	TFile *f_sim = NULL;
	if(TFile::Open(TString::Format("%s__exp.root", fname)) != NULL) f_exp = TFile::Open(TString::Format("%s__exp.root", fname));
	if(TFile::Open(TString::Format("simdir/%s__sim.root", fname)) != NULL) f_sim = TFile::Open(TString::Format("simdir/%s__sim.root", fname));
	if (_type.EqualTo("a") && (f_exp==NULL || f_sim==NULL)) {
		printf("f_exp AND f_sim needed for acceptance study\n");
		return;
	}else if (_type.EqualTo("d") && f_sim==NULL){
		printf("f_sim need for dead-cells study\n");
		return;
	}else if (!_type.EqualTo("a") && !_type.EqualTo("d")){
		printf("%s _type is not valid!\n",_type.Data());
		return;
	}
		
	//for sim stats
	ofstream f_simstats; 
	//f_simstats.open("simstats.txt");
		  
	/* *** from fname, determine number of nQ2bins, nQ2Wbins *** */
	TString fName = f_sim->GetName();
	TObjArray* fNameArray = fName.Tokenize("__");
	TString Q2bng = fNameArray->At(1)->GetName();
	TString Wbng  = fNameArray->At(2)->GetName();
	TString dtype = fNameArray->At(3)->GetName();
	TObjArray* Q2bngArray = Q2bng.Tokenize("-");
	TObjArray* WbngArray  = Wbng.Tokenize("-");
	Int_t nQ2bins = atoi(Q2bngArray->At(0)->GetName());
	Int_t nWbins = atoi(WbngArray->At(0)->GetName());
	Int_t nQ2Wbins = nQ2bins*nWbins;
	printf("nQ2bins:nWbins = %d:%d\n", nQ2bins, nWbins);
  
	/* *** Create Canvas and divide it according to nQ2bins, nWbins *** */
	TRandom *R = new TRandom(time(0));
	rNum = R->Integer(1000);
	TCanvas *c_accStudy[100];
	TCanvas *c_dcStudy;
		
	if(_type=="a"){
		for (Int_t iQ2Wbin=0;iQ2Wbin<nQ2Wbins;iQ2Wbin++){
			if (mDraw){
				TString name = TString::Format("%s_%d", "c_accStudy", iQ2Wbin+1);
				TString title = TString::Format("%s_%d","acceptance study", iQ2Wbin+1);
				c_accStudy[iQ2Wbin]= new TCanvas(name, title, 400, 300);
				c_accStudy[iQ2Wbin]->Divide(4,2);
			}
		}
	}else if(_type=="d"){
		if (mDraw) {
			TString name = TString::Format("%s", "c_dcStudy");
			TString title = TString::Format("%s","dead cells study");
			c_dcStudy = new TCanvas(name, title, 400, 300);
			c_dcStudy->Divide(nQ2Wbins);
			c_dcStudy->AddExec("tnail_dcStudy", "tnail_dcStudy(c_dcStudy, c_tnail_dcStudy)");
			TCanvas* c_tnail_dcStudy = new TCanvas("c_tnail_dcStudy", "tnail_dcStudy");
		}
		
		f_simstats.open(TString::Format("simstats_%s__sim.root", fname));
		
	}
		
	f_sim->GetListOfKeys()->Print(); //all Q2W dirs
	TIter nextkey(f_sim->GetListOfKeys());
	TKey *key;

	//Loop over Q2W dirs, get h5Ds and "linearize" them to h1Ds
	Int_t iQ2Wbin = 0;
	while (key = (TKey*)nextkey()) {
		Q2Wdirname = key->GetName();
		if(Q2Wdirname.Contains("hYW"))continue;
		printf("Q2Wdir = %s\n", Q2Wdirname.Data());

		h5Dpath[TH]        = TString::Format("%s/hY5D/Varset1/hY5D_TH", Q2Wdirname.Data());
		h5Dpath[RECO]      = TString::Format("%s/hY5D/Varset1/hY5D_RECO", Q2Wdirname.Data());
		h5Dpath[ACC]       = TString::Format("%s/hY5D/Varset1/hY5D_ACC", Q2Wdirname.Data());
		h5Dpath[ACC_CORR]  = TString::Format("%s/hY5D/Varset1/hY5D_ACC_CORR", Q2Wdirname.Data());
		h5Dpath[HOLE]  = TString::Format("%s/hY5D/Varset1/hY5D_HOLE", Q2Wdirname.Data());
		h5Dpath[FULL]  = TString::Format("%s/hY5D/Varset1/hY5D_FULL", Q2Wdirname.Data());
				
		//!Get h5Ds hists from file
		for (Int_t iSeq=0;iSeq<nSeq;iSeq++){
			seq_t seq = static_cast<seq_t>(iSeq);
			if (_type.EqualTo("a") && !(seq==TH) && !(seq==ACC)){
				 h5D[EXP][seq] = (THnSparse*)f_exp->Get(h5Dpath[seq].Data());
				 nBins_h5D[EXP][seq] = h5D[EXP][seq]->GetNbins();
			} 
			h5D[SIM][seq] = (THnSparse*)f_sim->Get(h5Dpath[seq].Data());
			nBins_h5D[SIM][seq] = h5D[SIM][seq]->GetNbins();
		}
		if      (_type.EqualTo("a")) nFBins = nBins_h5D[EXP][RECO];
		else if (_type.EqualTo("d")) nFBins = nBins_h5D[SIM][TH];
                       
		//!Create h1Ds
		makeHists_h1D(); 
			
		//! Linearize h5D --> h1D
		for (Int_t iFBin=0;iFBin<nFBins;iFBin++){
			Int_t bincoord[5] = {0};
			if (_type=="a"){
				//!first do EXP.RECO (obtain bincoord[])
				FillLinHist(h1D[EXP][RECO],     h5D[EXP][RECO],     iFBin, bincoord);
				FillLinHist(h1D[EXP][ACC_CORR], h5D[EXP][ACC_CORR], iFBin, bincoord); 
				FillLinHist(h1D[SIM][TH],       h5D[SIM][TH],       iFBin, bincoord);
				FillLinHist(h1D[SIM][RECO],     h5D[SIM][RECO],     iFBin, bincoord);
				FillLinHist(h1D[SIM][ACC],      h5D[SIM][ACC],      iFBin, bincoord);
				FillLinHist(h1D[SIM][ACC_CORR], h5D[SIM][ACC_CORR], iFBin, bincoord);
			}else if (_type=="d"){
				//!first do SIM.TH (obtain bincoord[])
				FillLinHist(h1D[SIM][TH],   h5D[SIM][TH],   iFBin, bincoord);
				FillLinHist(h1D[SIM][RECO], h5D[SIM][RECO], iFBin, bincoord);
			}
		}
		
		if(_type=="a"){
			for (Int_t iSeq=0;iSeq<nSeq;iSeq++){
				seq_t seq = static_cast<seq_t>(iSeq);
				if ( seq==HOLE || seq==FULL) continue;
				
				if( (seq!=TH) && (seq!=ACC) ){
					if (mDraw){
						c_accStudy[iQ2Wbin]->cd( (seq+1) );
						h1D[EXP][seq]->Draw("P");
					}
				}
				if (mDraw){
					c_accStudy[iQ2Wbin]->cd( (seq+1) + (4) );
					h1D[SIM][seq]->Draw("P");
				}
			}
		}else if(_type=="d"){		
			if (mDraw){
				c_dcStudy->cd(iQ2Wbin+1);
				h1D[SIM][TH]->Draw("e1p");
				h1D[SIM][RECO]->Draw("e1p sames");
			}
						
			Int_t n0simRECObins = GetNum0Bins(h1D[SIM][RECO]);
			//Int_t nFsimRECObins = h1D[SIM][RECO]->GetEntries() - n0simRECObins;
			Int_t nFsimRECObins = h1D[SIM][RECO]->GetNbinsX() - n0simRECObins;
			Int_t nRECOevts = h1D[SIM][RECO]->Integral();
			Int_t nTHevts = h1D[SIM][TH]->Integral();
			f_simstats << Q2Wdirname << "	" << n0simRECObins << "	" << nFsimRECObins << "	" << nRECOevts << "	" << nTHevts << endl;
		}
		
		iQ2Wbin+=1;
	}
	f_simstats.close();
	
}

void tnail_dcStudy(TCanvas* c_dcStudy, TCanvas* c_tnail_dcStudy) {
	TPad *sel = (TPad*)gPad->GetSelectedPad();
	int px = gPad->GetEventX();
	int py = gPad->GetEventY();
	if (sel && sel != c_tnail_dcStudy && sel != c_dcStudy) {
		c_tnail_dcStudy->cd();
		TPad *newpad = (TPad*)sel->Clone();
		c_tnail_dcStudy->GetListOfPrimitives()->Add(newpad);
		newpad->SetPad(0,0,1,1);
		//selold_tnail = newpad;
		c_tnail_dcStudy->Update();
		c_dcStudy->cd();
	}
}

void makeHists_h1D(){
	for (Int_t iSeq=0;iSeq<nSeq;iSeq++){
		seq_t seq = static_cast<seq_t>(iSeq);
		//!EXP
		if (!(seq == TH) && !(seq==ACC)){
			TString name = TString::Format("%s.%s_nEvtsVbin", "EXP", seqTitle[seq].Data());
			TString title = TString::Format("nEvts V bin: %s %s", seqTitle[seq].Data(), Q2Wdirname.Data());
			h1D[EXP][seq] = new TH1D(name, title, nFBins, 0.5, nFBins+0.5);
			h1D[EXP][seq]->SetMinimum(-10000);
			h1D[EXP][seq]->SetLineColor(kBlack);
			h1D[EXP][seq]->SetMarkerStyle(33);
			h1D[EXP][seq]->SetMarkerSize(0.5);
			h1D[EXP][seq]->SetMarkerColor(kBlue-7);
		}
		TString name = TString::Format("%s.%s_nEvtsVbin", "SIM", seqTitle[seq].Data());
		TString title = TString::Format("nEvts V bin: %s %s", seqTitle[seq].Data(), Q2Wdirname.Data());
		h1D[SIM][seq] = new TH1D(name, title, nFBins, 0.5, nFBins+0.5);
		h1D[SIM][seq]->SetMinimum(-10000);
		if(seq==ACC) {
			h1D[SIM][seq]->SetMinimum(0);
			h1D[SIM][seq]->SetMaximum(0.4);
		}
		if(_type=="a"){
			h1D[SIM][seq]->SetLineColor(kBlack);
			h1D[SIM][seq]->SetMarkerStyle(33);
			h1D[SIM][seq]->SetMarkerSize(0.5);
			if (seq==TH)                         h1D[SIM][seq]->SetMarkerColor(kGreen-7);
			else if (seq==RECO || seq==ACC_CORR) h1D[SIM][seq]->SetMarkerColor(kRed-7);
			else if (seq==ACC)                   h1D[SIM][seq]->SetMarkerColor(kOrange-7);
		}else if (_type=="d"){
			if (seq==TH)                         h1D[SIM][seq]->SetLineColor(kGreen-7);
			else if (seq==RECO || seq==ACC_CORR) h1D[SIM][seq]->SetLineColor(kRed-7);
			else if (seq==ACC)                   h1D[SIM][seq]->SetLineColor(kOrange-7);
		}
			
	}
}


void FillLinHist(TH1* h1, THnSparse* hN, Int_t iFBin, Int_t* FBincoord){
	Float_t binc;
	if (FBincoord[0]==0 && FBincoord[1]==0 && FBincoord[2]==0 && FBincoord[3]==0 && FBincoord[4]==0){
		binc   = hN->GetBinContent(iFBin, FBincoord);
		//Float_t binerr = hN->GetBinError(iFBin);
		h1->SetBinContent(iFBin+1, binc);
	}else{
		Int_t iBin   = hN->GetBin(FBincoord);
		binc = hN->GetBinContent(iBin);
		//Float_t binerr = hN->GetBinError(iBin);
		h1->SetBinContent(iFBin+1, binc);
	}
	if(_type=="a"){
		if      (binc!=0) h1->SetBinError  (iFBin+1, 0);
		else if (binc==0) h1->SetBinError  (iFBin+1, 10000);
	}else if(_type=="d"){
		if      (binc!=0) h1->SetBinError  (iFBin+1, TMath::Sqrt(binc));
		else if (binc==0) h1->SetBinError  (iFBin+1, 10000);
	}
}

Int_t GetNum0Bins(TH1* h1){
	Int_t n0Bins = 0;
	Int_t nBins = h1->GetNbinsX();
	for(Int_t iBin=0;iBin<nBins;iBin++){
		if (h1->GetBinContent(iBin+1)==0) n0Bins+=1;
	}
	return n0Bins;
}



