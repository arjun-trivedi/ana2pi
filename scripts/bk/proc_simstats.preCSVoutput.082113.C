#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TH1.h>

#include <iostream>

#include "proc_simstats.h"
#include "myTHnTool.h"

void proc_simstats(TString sim, TString xsectype="at"){
	myTHnTool hNtool(kFALSE); //set myTHnTool._dbgMode = kFALSE; no couts from myTHnTool::GetBinStats()
	Int_t nQ2Wbins; 
	Int_t nSims;
	TString q2wbng;
	TString* simDir;
	TString expDir;
	TString topNames[nTops];
	if(sim=="e1fs1"){
		nQ2Wbins = nQ2Wbins_e1fs1;
		nSims=nSims_e1fs1;
		q2wbng = q2wbng_e1fs1;
		simDir = simDir_e1fs1;
		expDir = expDir_e1fs1;
	}else if (sim=="e1fs2"){
		nQ2Wbins = nQ2Wbins_e1fs2;
		nSims=nSims_e1fs2;
		q2wbng = q2wbng_e1fs2;
		simDir = simDir_e1fs2;
		expDir = expDir_e1fs2;
	}else{
		printf("simulation could not be determined\n");
		return;
	}
	if(xsectype=="at"){
		memcpy(topNames, at_topNames, sizeof (at_topNames));
	}else if(xsectype=="vm"){
		memcpy(topNames, vm_topNames, sizeof (vm_topNames));
	}else{
		printf("xsectype could not be determined\n");
		return;
	}
	//! create h4[nTops][nSims][nQ2Wbins][nVarsets]
	TFile* fout = new TFile(TString::Format("simstats_%s.root", xsectype.Data()), "RECREATE");
	Int_t hdim = 4;
	Int_t nbins[4]      = {nTops,   nSims,   nQ2Wbins,   nVarsets};
	Double_t xmin[4]    = {1,       1,       1,          1};
	Double_t xmax[4]    = {nTops+1, nSims+1, nQ2Wbins+1, nVarsets+1};
	THnSparse* h4_n0sr  = new THnSparseF("h4_n0sr",  "h4_n0sr",  hdim, nbins, xmin, xmax);
	THnSparse* h4_n0acc = new THnSparseF("h4_n0acc", "h4_n0acc", hdim, nbins, xmin, xmax);
    THnSparse* h4_nFsr  = new THnSparseF("h4_nFsr",  "h4_nFsr",  hdim, nbins, xmin, xmax);
    THnSparse* h4_nFth  = new THnSparseF("h4_nFth", "h4_nFth", hdim, nbins, xmin, xmax);
	
	for (Int_t iTop=0;iTop<nTops;iTop++){//begin Top loop
		TFile* fexp = NULL;  
		fexp = new TFile( TString::Format("%s/%s__%s__exp.root", expDir.Data(), topNames[iTop].Data(), q2wbng.Data()) );
		if (fexp==NULL) {
			printf("%s does not exist\n", fexp->GetName());
		continue;
		}
                           
		for (Int_t iSim=0;iSim<nSims;iSim++){//begin Sim loop
			TFile* fsim = NULL;
			fsim = TFile::Open( TString::Format("%s/%s__%s__sim.root", simDir[iSim].Data(), topNames[iTop].Data(), q2wbng.Data()) );
			if (fsim==NULL) {
				printf("%s does not exist\n", fsim->GetName());
				continue;
			}
			//!Loop over Q2W dirs, get h5Ds and their yields
			Int_t iQ2Wbin=0;
			TIter nextkey(fsim->GetListOfKeys());
			TKey *key;
			while (key = (TKey*)nextkey()) {//begin Q2W loop
				TString Q2Wdirname = key->GetName();
				if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
				for(Int_t iVarset=0;iVarset<nVarsets;iVarset++){//begin Varset loop
					Double_t coord[4] = {iTop+1, iSim+1, iQ2Wbin+1, iVarset+1};
					char hname[100];
					sprintf(hname, "%s/hY5D/Varset%d/hY5D_TH", Q2Wdirname.Data(),iVarset+1);
					THnSparse* hY5D_TH = (THnSparse*)fsim->Get(hname);
					sprintf(hname, "%s/hY5D/Varset%d/hY5D_RECO", Q2Wdirname.Data(),iVarset+1);
					THnSparse* hY5D_simRECO = (THnSparse*)fsim->Get(hname);
					THnSparse* hY5D_expRECO = (THnSparse*)fexp->Get(hname);
					sprintf(hname, "%s/hY5D/Varset%d/hY5D_ACC", Q2Wdirname.Data(),iVarset+1);
					THnSparse* hY5D_ACC = (THnSparse*)fsim->Get(hname);

					h4_n0acc->Fill(coord, hNtool.GetNbinsEq0(hY5D_ACC,hY5D_expRECO));
					h4_n0sr ->Fill(coord, hNtool.GetNbinsEq0(hY5D_simRECO,hY5D_TH));
					h4_n0acc->GetAxis(2)->SetBinLabel(iQ2Wbin+1, Q2Wdirname.Data()); //bin labels for Q2W bins
					h4_n0sr ->GetAxis(2)->SetBinLabel(iQ2Wbin+1, Q2Wdirname.Data());

					h4_nFsr->Fill(coord, hNtool.GetNbinsNotEq0(hY5D_simRECO));
					h4_nFth->Fill(coord, hNtool.GetNbinsNotEq0(hY5D_TH));
					h4_nFsr->GetAxis(2)->SetBinLabel(iQ2Wbin+1, Q2Wdirname.Data()); //bin labels for Q2W bins
					h4_nFth->GetAxis(2)->SetBinLabel(iQ2Wbin+1, Q2Wdirname.Data());
				}//end Varset loop
			iQ2Wbin+=1;	
			}//end Q2W bin loop
		}//end Sim loop
	}//end Top loop
	fout->cd();
	
    h4_n0acc->Write();
	h4_n0sr->Write();
        
	h4_nFsr->Write();
	h4_nFth->Write(); 

	fout->Write();
	fout->Close();
}
