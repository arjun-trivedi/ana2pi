/*
[11-22-13]
In a CSV format, store the following tabular data:
"Sim","Top","Q2Wbin","Varset","nFB_ST","nEB_SR","nFB_SR","nEB_SA"

where:
nFB_ST = Total number of Bins that are Filled (FB) by `genev` (ST) in the 2pi Reaction Phase Space (PS)
nEB_SR = Out of nFB_ST, the number of Bins that DO NOT HAVE Reconstructed events.
nFB_SR = Out of nFB_ST, the number of Bins that HAVE Reconstructed events.
nEB_SA = Out of nFB_ER (number of Bins that are Filled by Reconstructed events in Experiment), the number
         of bins that have no Acceptance.

Hypothesis for "Complete" Simulation
====================================
Before stating the hypothesis, I want to note the following 2 important assumptions:
1. It is assumed that `genev`=Nature and therefore, `genev` throws events within the Natural PS
2. It is assumed that `GSIM`= CLAS Detector

Counting nEB_SR/nFB_SR is the primary means for keeping track of when the Simulation is "complete". Ideally,
for all nFB_ST that lie within the Fiducial Volume of the CLAS detector, there should be a FB_SR. Any of 
the nFB_ST events that lie outside this Fiducial volume, should count for nEB_SR.
Therefore:
nFB_SR = nFB_ST-nEB-SR

In addition, nEB_SA can serve as another (independent?) gauge. For if the Simulation is "Complete",
nEB_SA should approach its asymptotic limit too. Ideally, if ER data has no Background and the Simulation is 
ideal, this asymptotic limit should approach 0. However, the limit will approach a number that may be representative
of the number of bins in ER data that is filled with Background events.

Hence, once Simulation is complete:
1. nFB_ST,nFB_SR,nEB_SR,nEB_SA should approach their asymptotic limits
2. nFB_SR = nFB_ST-nEB-SR
*/
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TH1.h>
#include <TTimeStamp.h>

#include <iostream>
#include <fstream>
using namespace std;

#include "proc_simstats.h"
#include "myTHnTool.h"

void proc_simstats(TString sim, hel_t hel=UNPOL,TString xsectype="vm"){
	myTHnTool hNtool(kFALSE); //set myTHnTool._dbgMode = kFALSE; no couts from myTHnTool::GetBinStats()
	Int_t nq2wbins; 
	Int_t nsims;
	TString q2wbng;
	TString* simdirs;
	TString expdir;
	TString topnames[nTOP+1];
	if(sim=="e1fs1"){
		nsims    = e1fs1_nsims;
		simdirs  = e1fs1_simdirs;
		expdir   = e1fs1_expdir;
		nq2wbins = e1fs1_nq2wbins;
		q2wbng   = e1fs1_q2wbng;
	}else if (sim=="e1fs2"){
		nsims    = e1fs2_nsims;
		simdirs  = e1fs2_simdirs;
		expdir   = e1fs2_expdir;
		nq2wbins = e1fs2_nq2wbins;
		q2wbng   = e1fs2_q2wbng;
	}else{
		printf("simulation could not be determined\n");
		return;
	}
	if(xsectype=="at"){
		memcpy(topnames, at_topnames, sizeof (at_topnames));
	}else if(xsectype=="vm"){
		memcpy(topnames, vm_topnames, sizeof (vm_topnames));
	}else{
		printf("xsectype could not be determined\n");
		return;
	}
	
	ofstream fout;
	/*if      (hel==UNPOL)fout.open(TString::Format("simstats_%s.csv", xsectype.Data()));
	else if (hel==POS)  fout.open(TString::Format("simstats_%s_POS.csv", xsectype.Data())); 
	else if (hel==NEG)  fout.open(TString::Format("simstats_%s_NEG.csv", xsectype.Data()));*/
	//char fout_name[100];
	TString fout_name;
	if      (hel==UNPOL) fout_name = TString::Format("%s/simstats_%s.csv", expdir.Data(),xsectype.Data());
	else if (hel==POS)   fout_name = TString::Format("%s/simstats_%s_POS.csv", expdir.Data(),xsectype.Data()); 
	else if (hel==NEG)   fout_name = TString::Format("%s/simstats_%s_NEG.csv", expdir.Data(), xsectype.Data());
	FileStat_t fstat;
	if (gSystem->GetPathInfo(fout_name,fstat)!=0) {//i.e. fout_name does not exist
		fout.open(fout_name);
	}else{ //if fout_name exists then back it up and then create new fout_name
		TTimeStamp ts;
		UInt_t tstmp[6]; //year, month, day, hr,min, sec
		bool inUTC=kFALSE;
		int secOffset = 0;
		ts.GetDate(inUTC,secOffset,&tstmp[0],&tstmp[1],&tstmp[2]);
		ts.GetTime(inUTC,secOffset,&tstmp[3],&tstmp[4],&tstmp[5]);
		TString str_tstmp = TString::Format("%04d-%02d-%02d_%02d-%02d", tstmp[0],tstmp[1],tstmp[2],tstmp[3],tstmp[4]); 
		cout << str_tstmp << endl;
		TString fout_name_bk = TString::Format("%s_bk_%s", fout_name.Data(), str_tstmp.Data());
		TString shellcmd = TString::Format("cp -p %s %s", fout_name.Data(), fout_name_bk.Data());
		gSystem->Exec(shellcmd);
		fout.open(fout_name);
	}
	
	char fout_col_names[100];
	sprintf(fout_col_names,"%s,%s,%s,%s,%s,%s,%s,%s\n","Sim","Top","Q2Wbin","Varset","nFB_ST","nEB_SR","nFB_SR","nEB_SA");
	fout << fout_col_names;
	
	for (Int_t iSim=0;iSim<nsims;iSim++){//begin Sim loop
		for (Int_t iTop=0;iTop<nTOP+1;iTop++){//begin Top loop
			TFile* fsim = NULL;
			TFile* fexp = NULL;
			fsim = TFile::Open( TString::Format("%s/%s__%s__sim.root", simdirs[iSim].Data(), topnames[iTop].Data(), q2wbng.Data()) );
			if     (hel==UNPOL)fexp = TFile::Open( TString::Format("%s/%s__%s__exp.root", expdir.Data(),        topnames[iTop].Data(), q2wbng.Data()) );
			else if(hel==POS||hel==NEG)  fexp = TFile::Open( TString::Format("%s/%s__%s__pol__exp.root", expdir.Data(),        topnames[iTop].Data(), q2wbng.Data()) ); 
			if (fsim==NULL||fexp==NULL) {
				if (fsim==NULL)printf("%s does not exist\n", fsim->GetName());
				if (fexp==NULL)printf("%s does not exist\n", fsim->GetName());
				continue;
			}

			//!Loop over Q2W dirs, get h5Ds and their yields
			Int_t iQ2Wbin=0;
			TIter nextkey(fsim->GetListOfKeys());
			TKey *key;
			while (key = (TKey*)nextkey()) {//begin Q2W loop
				TString Q2Wdirname = key->GetName();
				if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
				for(Int_t iVarset=0;iVarset<nVARSET;iVarset++){//begin Varset loop
					printf("Obtaining simstats for Sim:Top:Q2Wbin:Varset %d:%d:%s:%d \n", iSim+1,iTop+1,Q2Wdirname.Data(),iVarset+1);
					char hname[100];

					sprintf(hname, "%s/hY5D/Varset%d/hY5D_TH", Q2Wdirname.Data(),iVarset+1);
					THnSparse* hY5D_TH = (THnSparse*)fsim->Get(hname);

					sprintf(hname, "%s/hY5D/Varset%d/hY5D_RECO", Q2Wdirname.Data(),iVarset+1);
					THnSparse* hY5D_SR = (THnSparse*)fsim->Get(hname);

					THnSparse* hY5D_ER = NULL;
					if(hel==UNPOL)    hY5D_ER = (THnSparse*)fexp->Get(hname);
					else if(hel==POS) hY5D_ER = (THnSparse*)fexp->Get(TString::Format("%s/hY5D_POS/Varset%d/hY5D_RECO", Q2Wdirname.Data(),iVarset+1));
					else if(hel==NEG) hY5D_ER = (THnSparse*)fexp->Get(TString::Format("%s/hY5D_NEG/Varset%d/hY5D_RECO", Q2Wdirname.Data(),iVarset+1));

					sprintf(hname, "%s/hY5D/Varset%d/hY5D_ACC", Q2Wdirname.Data(),iVarset+1);
					THnSparse* hY5D_ACC = (THnSparse*)fsim->Get(hname);

					int nFB_ST = hNtool.GetNbinsNotEq0(hY5D_TH);
					int nEB_SR = hNtool.GetNbinsEq0(hY5D_SR,hY5D_TH);
					int nFB_SR = hNtool.GetNbinsNotEq0(hY5D_SR,hY5D_TH);
					int nEB_SA = hNtool.GetNbinsEq0(hY5D_ACC,hY5D_ER);
					
					char fout_data[100];
					sprintf(fout_data,"%d,%d,%d,%d,%d,%d,%d,%d\n",iSim+1,iTop+1,iQ2Wbin+1,iVarset+1,nFB_ST,nEB_SR,nFB_SR,nEB_SA);
					fout << fout_data;
				}//end Varset loop
			iQ2Wbin+=1;	
			}//end Q2W bin loop
		}//end Top loop
	}//end Sim loop

	printf("Closing output file...\n");
	fout.close();
	printf("Done Closing output file\n");
	return;
}
