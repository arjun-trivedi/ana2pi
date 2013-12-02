#define SelH10_cxx

#include "sel_h10.h"
#include <TH2.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TCollection.h>
#include <TObjString.h>
#include <TRegexp.h>
#include <TSystem.h>
#include "proc_eid.h"
#include "proc_efid.h"
#include "proc_pid.h"
#include "proc_skim_q.h"
#include "proc_mom_cor.h"
#include "proc_top.h"
#include "proc_skim_q2w.h"
#include "proc_fill_skim.h"
#include "proc_copy_h10.h"
#include "particle_constants.h"
#include "data_h10.h"

using namespace ParticleConstants;

SelH10::SelH10(TTree * /*tree*/) : TSelector()
{
	Info("SelH10", "");
	
	fFileOut = NULL;
	fProofFile = NULL;
	fOutFileName = NULL;
	
	Info("SelH10", "Done\n");
}

void SelH10::SetInputList(TList *input) {
	Info("SetInputList", "");
	fInput = input;
	if(fInput){
		Info("SetInputList", "fInput for TSelector set to %s", fInput->First()->GetTitle());
	}else{
		Info("SetInputList", "fInput not set!. Exiting\n");
		return;
	}
	if (fInput->FindObject("OUTPUTFILE"))       Info("SetInputList", "OUTPUTFILE set\n");
    if (fInput->FindObject("PROOF_OUTPUTFILE")) Info("SetInputList", "PROOF_OUTPUTFILE set\n");
    if (fInput->FindObject("H10TYPE")) Info("SetInputList", "H10TYPE set\n");
    
	Info("SetInputList", "Done\n");
	return;
}

void SelH10::Begin(TTree* /*tree*/)
{
	Info("Begin", "");
	Info("Begin", "Done...\n");
}

void SelH10::SlaveBegin(TTree* /*tree*/) 
{
	Info("SlaveBegin", "");
	
	//! Determine output file name; OUTPUTFILE/PROOF_OUTPUTFILE !//   
    TNamed *outFile      = (TNamed *) fInput->FindObject("OUTPUTFILE");
    TNamed *outFileProof = (TNamed *) fInput->FindObject("PROOF_OUTPUTFILE");
    if(outFile){
		fFileOut = new TFile(outFile->GetTitle(), "RECREATE");
		if (!fFileOut) {
			Info("SlaveBegin", "could not create '%s': instance is invalid!", fFileOut->GetName());
			Info("SlaveBegin", "Exiting\n");
			return;
		}
	}else if(outFileProof){// The file for merging
		fProofFile = new TProofOutputFile(outFileProof->GetTitle(), "M");
		fProofFile->SetOutputFileName(outFileProof->GetTitle());
		TDirectory *savedir = gDirectory;
		fFileOut = fProofFile->OpenFile("RECREATE");
		if (fFileOut && fFileOut->IsZombie()) SafeDelete(fFileOut);
		savedir->cd();

		// Cannot continue
		if (!fFileOut) {
			Info("SlaveBegin", "could not create '%s': instance is invalid!", fProofFile->GetName());
			Info("SlaveBegin", "Exiting\n");
			return;
		}
	}else{
		Info("SlaveBegin", "Output file not specified!");
		Info("SlaveBegin", "Exiting\n");
		return;
	}
	Info("SlaveBegin", "output in file '%s'", fFileOut->GetName());
	/* *** */
	
	
	//treenumber = -1; //atrivedi: --> Init() for Prooflite compatibility
	firstFileOfRun = kTRUE;
	fFileOut->cd();
	/* *** */
	
	//! Determine h10type & setup DataH10, DataAna & Procs !//
	h10type = (TNamed *) fInput->FindObject("H10TYPE");
	if(h10type){
		TObjArray *h10type_tokens = ((TString)h10type->GetTitle()).Tokenize(":");
		if(h10type_tokens->GetEntries() >= 2){
			// Determint h10type
			is_h10e1f = is_h10e16 = is_h10exp = is_h10sim = kFALSE;
			TObjArray *h10type_tokens = ((TString)h10type->GetTitle()).Tokenize(":");
			h10exp   = h10type_tokens->At(0)->GetName();
			h10dtype = h10type_tokens->At(1)->GetName();
			if (h10type_tokens->GetEntries() == 3) h10skim = h10type_tokens->At(2)->GetName();
			
			if      (h10exp.EqualTo("e1f")) is_h10e1f = kTRUE;
			else if (h10exp.EqualTo("e16")) is_h10e16 = kTRUE;
			else {
				Info("SlaveBegin", "Could not determine h10type.experiment!\n");
				Info("SlaveBegin", "!! CODE IS NOT DESIGNED TO WORK FROM THIS POINT ONWARDS !!");
			}
	
			if (h10dtype.EqualTo("exp"))      is_h10exp = kTRUE;
			else if (h10dtype.EqualTo("sim")) is_h10sim = kTRUE;
			else {
				Info("SlaveBegin", "Could not determine h10type.dtype!\n");
				Info("SlaveBegin", "!! CODE IS NOT DESIGNED TO WORK FROM THIS POINT ONWARDS !!");
			}
	
			Info("SlaveBegin","Going to set up DataH10 and Proc<X> with h10type: %s:%s:%s\n", h10exp.Data(), h10dtype.Data(), h10skim.Data());
			
			//Setup DataH10
			dH10 = new DataH10(h10type->GetTitle());
			
			//Setup DataAna(passed to each Proc) & Procs
			dAna = new DataAna(); 
			SetupProcs();
			
		}else{
			Info("SlaveBegin", "h10type identifier does not have enough information! Tokens < 2\n");
			Info("SlaveBegin", "!! CODE IS NOT DESIGNED TO WORK FROM THIS POINT ONWARDS !!");
		}
	}else {
		Info("SlaveBegin", "H10TYPE not set in InputList of SelH10. Going to call SelH10::Terminate()");
		Info("SlaveBegin", "!! CODE IS NOT DESIGNED TO WORK FROM THIS POINT ONWARDS !!");
		Terminate();
	}
	    
	/* *** */
	eventnum = numfilled = 0;
	blocksize = 100000;
	swmain.Start();
	swgroup.Start();
	
	Info("SlaveBegin", "Done\n");
	return;
	
}

void SelH10::Init(TTree *tree) {
	Info("Init", "");
		
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set branch addresses and branch pointers
	if (!tree) {
		Info("Init", "No TChain to process! Exiting\n");
		return;
	}
	fChain = tree;
	fChain->SetMakeClass(1);
	
	Info("Init", "Getting number of entries in TChain...");
	numInChain += fChain->GetEntries();
	Info("Init", "Number of entries in TChain = %d", numInChain);
	
	/* *** atrivedi 112912: Activate only used branches*** */
	Info("Init", "lastProc = %s", lastProcName->GetName());
	Info("Init", "Going to activate only selected branches");
	fChain->SetBranchStatus("*",0); //disable all branches
	
	//Basic: Branches needed for Basic running of ana2pi (all implemented Proc<Cuts/Corrections> and ProcTop
	fChain->SetBranchStatus("evthel", 1); //e1f
	fChain->SetBranchStatus("gpart", 1);
	fChain->SetBranchStatus("id", 1);
	fChain->SetBranchStatus("stat", 1);
	fChain->SetBranchStatus("p", 1);
	fChain->SetBranchStatus("cx", 1);
	fChain->SetBranchStatus("cy", 1);
	fChain->SetBranchStatus("cz", 1);
	fChain->SetBranchStatus("q", 1);
	fChain->SetBranchStatus("dc", 1);
	fChain->SetBranchStatus("cc", 1);
	fChain->SetBranchStatus("sc", 1);
	fChain->SetBranchStatus("ec", 1);
	fChain->SetBranchStatus("dc_part", 1); //atrivedi 042413
	fChain->SetBranchStatus("dc_stat", 1);
	fChain->SetBranchStatus("sc_part", 1); //atrivedi 042413
	fChain->SetBranchStatus("sc_sect", 1);
	fChain->SetBranchStatus("sc_pd", 1);
	fChain->SetBranchStatus("ec_part", 1); //atrivedi 042413
	fChain->SetBranchStatus("etot", 1);
	if (is_h10sim){
		fChain->SetBranchStatus("mcnentr",1);
		fChain->SetBranchStatus("mcid",1);
		fChain->SetBranchStatus("mcp",1);
		fChain->SetBranchStatus("mctheta",1);
		fChain->SetBranchStatus("mcphi",1);
		fChain->SetBranchStatus("mcm",1);
	}
		
	//Non-Basic: Additional branches needed when running ana2pi Monitoring mode
	fChain->SetBranchStatus("sc_r", 1);
	fChain->SetBranchStatus("sc_t", 1);
		
	//added for proc_eidmon (using at-eid)
	fChain->SetBranchStatus("dc_xsc", 1);
	fChain->SetBranchStatus("dc_ysc", 1);
	fChain->SetBranchStatus("dc_zsc", 1);
	fChain->SetBranchStatus("dc_cxsc", 1);
	fChain->SetBranchStatus("dc_cysc", 1);
	fChain->SetBranchStatus("dc_czsc", 1);
	fChain->SetBranchStatus("cc_segm", 1);
	fChain->SetBranchStatus("nphe", 1);
	fChain->SetBranchStatus("ec_eo", 1);
	fChain->SetBranchStatus("ec_ei", 1);
		
	//added for procfid_mon
	fChain->SetBranchStatus("ech_x", 1);
	fChain->SetBranchStatus("ech_y", 1);
	fChain->SetBranchStatus("ech_z", 1);
	/* *** */
	
	treenumber = -1;

	dH10->Bind(fChain);
			
	Info("Init", "Done\n");
	return;
	
}

Bool_t SelH10::Notify() {
	Info("Notify", "");
		
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.
	//Info("Notify", "Done\n");
	return kTRUE;
}

Int_t SelH10::GetEntry(Long64_t entry, Int_t getall) {
	Int_t retval = 0;
	if (fChain) {
		retval = fChain->GetTree()->GetEntry(entry, getall);
		//cout <<"Numer of bytes read for entry " << entry << " = " << retval << endl;
		//Reconcile();
		if(is_h10e16 || is_h10sim) {
			//Info("Init", "isE1fSim; going to call Reconcile()\n");
			dH10->Reconcile();
		}
		//return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
	}
}

Bool_t SelH10::Process(Long64_t entry)
{
	//Info("Process", "");
	dH10->Clear();
	dAna->Clear();
	if ( ++eventnum % blocksize == 0 ) {
		Float_t gtime = swgroup.RealTime();
		Float_t ttime = swmain.RealTime();
		Float_t percentProcessed = (Float_t)eventnum/numInChain*100;
		Float_t remaining = (100/percentProcessed*ttime-ttime)/60;
		printf("(%.2f\%) %lld/%.2f = %i events/sec | block = %i events/sec ... %.1f min remaining\n",percentProcessed,eventnum,ttime,((Int_t)(eventnum/ttime)),(Int_t)(blocksize/gtime),remaining);
		swgroup.Start();
		swmain.Start(kFALSE);
	}
	
	//cout << "opened file " << fChain->GetCurrentFile()->GetName() << endl;
	GetEntry(entry);
	
	if (treenumber != fChain->GetTreeNumber()) {
		Int_t prevRun = dH10->run;
		TChain *c = (TChain*) fChain;
		treenumber = c->GetTreeNumber();
		TString fn = c->GetCurrentFile()->GetName(); // was earlier GetFile() which crashed when using Prooflite
		memcpy(dH10->fn,fn.Data(),256);
		Int_t filenumber = c->GetTreeNumber()+1; //filenumber++;
		TRegexp _regexp_run("[0-9][0-9][0-9][0-9][0-9]");
		TString rnstr = fn(_regexp_run);
		dH10->run = rnstr.Atoi();
		Info("Process", " opened file #%d: %s (rn #%d)\n",filenumber,fn.Data(),dH10->run);
		fFileOut->cd();
		firstFileOfRun = (prevRun != dH10->run);
	}
			
	vector<EpProcessor*>::iterator itrProcs;
	for (itrProcs = procs.begin(); itrProcs < procs.end(); itrProcs++){
		(*itrProcs)->pass = kFALSE;
	}
	topProc->handle(dH10);
	if (procs.back()->pass){
		numfilled++;
	}
	//Info("Process", "Done\n");
	return kTRUE;
}

void SelH10::SlaveTerminate()
{
	Info("SlaveTerminate", "");
	// The SlaveTerminate() function is called after all entries or objects
	// have been processed. When running with PROOF SlaveTerminate() is called
	// on each slave server.
	
	if (fProofFile) {
		Bool_t cleanup = kTRUE;
		TDirectory *savedir = gDirectory;
		
		fFileOut->cd();
		Info("SlaveTerminate","Going to write hists from procs");
		Info("SlaveTerminate","last processor = %s", lastProcName->GetName());
		vector<EpProcessor*>::iterator itrProcs;
		for (itrProcs = procs.begin(); itrProcs < procs.end(); itrProcs++){
			(*itrProcs)->write();
		}
		
		cleanup = kFALSE;
		gDirectory = savedir;
		fFileOut->Close();
		// Cleanup or register
		if (cleanup) {
			Info("SlaveTerminate", "nothing to save: just cleanup everything ...");
			TUrl uf(*(fFileOut->GetEndpointUrl()));
			SafeDelete(fFileOut);
			gSystem->Unlink(uf.GetFile());
			SafeDelete(fProofFile);
			delete dH10;
			delete dAna;
		} else {
			Info("SlaveTerminate", "objects saved into '%s%s': sending related TProofOutputFile ...",
					fProofFile->GetFileName(), fProofFile->GetOptionsAnchor());
			fProofFile->Print();
			fOutput->Add(fProofFile);
			delete dH10;
			delete dAna;
		}
	}else {
		Info("SlaveTerminate", "Going to write non-Proof file %s", fFileOut->GetName());
		fFileOut->Write();
		fFileOut->Close();
		delete dH10;
		delete dAna;
	}
	
	Float_t ttime = swmain.RealTime();
	Float_t percentProcessed = (Float_t)eventnum/numInChain*100;
	Info("SlaveTerminate","Total: (%.2f\%) %lld/%.2f = %i events/sec",percentProcessed,eventnum,ttime,(Int_t)(eventnum/ttime));
	//while (!procs.empty()) procs.pop_back();
	//delete lastProcName;
			
	Info("SlaveTerminate", "Done\n");
	return;

}

void SelH10::Terminate()
{
	Info("Terminate", "");
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.
	
	Info("Terminate", "Done\n");
	return;
	
}

TDirectory* SelH10::mkdir(const char* dirname)
{
	TString dname = dirname;
	TDirectory *ret = 0;
	Int_t numTries = 0;
	while ( ret == 0 && numTries++<10 ) {
		if (numTries>1) dname += numTries;
		ret = fFileOut->mkdir(dname.Data());
	}
	return ret;
}

void SelH10::SetupProcs(){
	//Get list of processors
	TString option = GetOption();
	TCollection *procList = option.Tokenize(":");
	Info("SlaveBegin", "List of procs received:\n");
	procList->Print();
			
	//instantiate topProc; by design, topProc "builds" a cascaded chain of Procs
	topProc = new EpProcessor();
				
	if (procList->IsEmpty()) {
		Info("SlaveBegin", "no procs specified; building default pipeline\n");
		EpProcessor *proc;
		proc = new ProcEid(fFileOut->mkdir("eid"), dAna, h10type->GetTitle());
		procs.push_back(proc);
		topProc->add(proc);
		lastProcName= new TObjString("eid");
		Info("SlaveBegin","last processor = %s", lastProcName->GetName());
	} else {
		TIter iter(procList);
		while(TObjString *objStr = (TObjString*)iter.Next()) {
			EpProcessor *proc;
			TString str = objStr->GetString();
			if (str.EqualTo("eid")) 	        proc = new ProcEid(mkdir("eid"), dAna, h10type->GetTitle());
			else if (str.EqualTo("eidmon"))     proc = new ProcEid(mkdir("eid"), dAna, h10type->GetTitle(),kTRUE);
			else if (str.EqualTo("eidmononly")) proc = new ProcEid(mkdir("eid"), dAna, h10type->GetTitle(),kTRUE,kTRUE);
			else if (str.EqualTo("efid"))       proc = new ProcEFid(mkdir("fid"), dAna, h10type->GetTitle());
			else if (str.EqualTo("efidmon"))    proc = new ProcEFid(mkdir("fid"), dAna, h10type->GetTitle(),kTRUE);
			else if (str.EqualTo("efidmononly"))proc = new ProcEFid(mkdir("fid"), dAna, h10type->GetTitle(),kTRUE,kTRUE);
			else if (str.EqualTo("qskim")) 	    proc = new ProcSkimQ(mkdir("qskim"), dAna, h10type->GetTitle());
			else if (str.EqualTo("mom"))	    proc = new ProcMomCor(mkdir("mom"), dAna, h10type->GetTitle());
			else if (str.EqualTo("pid")) 	    proc = new ProcPid(mkdir("pid"), dAna, h10type->GetTitle());
			else if (str.EqualTo("pidmon"))     proc = new ProcPid(mkdir("pid"), dAna, h10type->GetTitle(),kTRUE);
			else if (str.EqualTo("pidmononly")) proc = new ProcPid(mkdir("pid"), dAna, h10type->GetTitle(),kTRUE,kTRUE);
			else if (str.EqualTo("top"))	    proc = new ProcTop(mkdir("top"), dAna, h10type->GetTitle());
			else if (str.EqualTo("q2wskim"))	proc = new ProcSkimQ2W(mkdir("q2wskim"), dAna, h10type->GetTitle());
			else if (str.EqualTo("fillskim"))	proc = new ProcFillSkim(mkdir("skim"), dAna, h10type->GetTitle());
			else if (str.EqualTo("copyh10"))	proc = new ProcCopyH10(fFileOut, dAna, h10type->GetTitle());
			else {
				Info("Init","%s unrecognized processor\n",str.Data());
				continue;
			}
			procs.push_back(proc);
			topProc->add(proc);
			lastProcName = objStr; //works because after last iteration, lastProcName is not updated
		}
		Info("SlaveBegin","last processor = %s", lastProcName->GetName());
	}
}
