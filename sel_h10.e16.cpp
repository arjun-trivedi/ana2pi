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
#include "proc_skim_q.h"
#include "proc_tops.h"
#include "proc_mom_cor.h"
#include "proc_fid.h"
#include "particle_constants.h"

using namespace ParticleConstants;

SelH10::SelH10(TTree * /*tree*/) : TSelector()
{
	Info("SelH10", "");
	
	tH10skim = NULL;
	
	fFileOut = NULL;
	fProofFile = NULL;
	fOutFileName = NULL;
	
	Info("SelH10", "Exiting\n");
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
    
	Info("SetInputList", "Exiting\n");
	return;
}

void SelH10::Begin(TTree * /*tree*/)
{
	Info("Begin", "");
	Info("Begin", "Exiting...\n");
}

void SelH10::SlaveBegin(TTree * /*tree*/)
{
	Info("SlaveBegin", "");
	   
    TString option = GetOption();
    
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
	
	//treenumber = -1; //atrivedi: --> Init() for Prooflite compatibility
	firstFileOfRun = kTRUE;
	fFileOut->cd();
	tOmega = new TTree("o","#omega topologies",1);
	tOmega->Branch("event",&data,16000,1);
	tQnorm = new TTree("fcq","luminosity block info per event");
	tQnorm->Branch("q_lb",&fLumBlock.q_lb,"q_lb/F");
	tQnorm->Branch("nEvnt",&fLumBlock.nEvnt,"nEvnt/F");
	tLumBlock = new TTree("lb","Luminosity Blocks",1);
	tLumBlock->Branch("run",&fLumBlock.run,"run/I");
    tLumBlock->Branch("event_first", &fLumBlock.evnt_first, "event_first/I");
    tLumBlock->Branch("event_last", &fLumBlock.evnt_last, "event_last/I");
    tLumBlock->Branch("id_last", &fLumBlock.id_last, "id_last/I");
    tLumBlock->Branch("nEvnt",&fLumBlock.nEvnt,"nEvnt/I");
    tLumBlock->Branch("nEvnt_e",&fLumBlock.nEvnt_e,"nEvnt_e/I");
    tLumBlock->Branch("nEvnt_ep",&fLumBlock.nEvnt_ep,"nEvnt_ep/I");
    tLumBlock->Branch("nEvnt_eppip",&fLumBlock.nEvnt_eppip,"nEvnt_eppip/I");
    tLumBlock->Branch("nEvnt_eppim",&fLumBlock.nEvnt_eppim,"nEvnt_eppim/I");
    tLumBlock->Branch("nEvnt_eppippim",&fLumBlock.nEvnt_eppippim,"nEvnt_eppippim/I");
    tLumBlock->Branch("q_lb",&fLumBlock.q_lb,"q_lb/F");
    tLumBlock->Branch("q_l",&fLumBlock.q_l,"q_l/F");
    tLumBlock->Branch("t_l",&fLumBlock.t_l,"t_l/F");
	TCollection *procs = option.Tokenize(":");
	if (procs->IsEmpty()) {
		Info("SlaveBegin", "no processors specified; building default pipeline\n");
		EpProcessor *proc;
		proc = new ProcEid(fFileOut->mkdir("eid"));
		processors.push_back(proc);
		top.add(proc);
		lastProcName= new TObjString("eid");
		Info("SlaveBegin","last processor = %s", lastProcName->GetName());
	} else {
		TIter iter(procs);
		while(TObjString *objStr = (TObjString*)iter.Next()) {
			EpProcessor *proc;
			TString str = objStr->GetString();
			if ( str.EqualTo("eid") ) proc = new ProcEid(mkdir("eid"));
			else if ( str.EqualTo("mom") ) proc = new ProcMomCor(mkdir("mom"));
			else if ( str.EqualTo("qskim") ) proc = new ProcSkimQ(mkdir("qskim"));
			//else if ( str.EqualTo("pid") ) proc = new ProcPid(fFileOut->mkdir("pid"));
			else if ( str.EqualTo("tops") ) proc = new ProcTops(mkdir("tops"));
			else if ( str.EqualTo("fid") ) proc = new ProcFid(mkdir("fid"));
			else {
				printf("%s unrecognized processor\n",str.Data());
				continue;
			}
			processors.push_back(proc);
			top.add(proc);
			lastProcName = objStr; //works because after last iteration, lastProcName is not updated
		}
		Info("SlaveBegin","last processor = %s", lastProcName->GetName());
	}
	eventnum = numfilled = 0;
	blocksize = 100000;
	swmain.Start();
	swgroup.Start();
	
	Info("SlaveBegin", "Exiting\n");
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
	
	//atrivedi
	if(!tH10skim){
		fFileOut->cd();
		tH10skim = (TTree*) fChain->GetTree()->CloneTree(0);
	}
	
	isSim = kTRUE;
	TBranch *testbranch = fChain->GetBranch("rf_time"); //e16
	if (testbranch) isSim = kFALSE;
	isExp = !isSim;
	data.h10.is_sim = isSim;
	
	/* *** atrivedi 112912: Activate only used branches*** */
	Info("Init", "lastProc = %s", lastProcName->GetName());
	if (strcmp(lastProcName->GetName(), "tops") == 0){ //atrivedi 012013: i.e. only when making yields
		Info("Init", "Going to activate only selected branches");
		fChain->SetBranchStatus("*",0); //disable all branches
		//main branches
		fChain->SetBranchStatus("id", 1);
		fChain->SetBranchStatus("p", 1);
		fChain->SetBranchStatus("cx", 1);
		fChain->SetBranchStatus("cy", 1);
		fChain->SetBranchStatus("cz", 1);
		fChain->SetBranchStatus("q", 1);
		fChain->SetBranchStatus("sc", 1);
		fChain->SetBranchStatus("dc", 1);
		fChain->SetBranchStatus("ec", 1);
		fChain->SetBranchStatus("cc", 1);
		fChain->SetBranchStatus("gpart", 1);
		fChain->SetBranchStatus("sc_sect", 1);
		fChain->SetBranchStatus("sc_pd", 1);
		//fChain->SetBranchStatus("evthel", 1); //e16
		//added for proc_eid (using at-eid)
		fChain->SetBranchStatus("dc_xsc", 1);
		fChain->SetBranchStatus("dc_ysc", 1);
		fChain->SetBranchStatus("dc_zsc", 1);
		fChain->SetBranchStatus("dc_cxsc", 1);
		fChain->SetBranchStatus("dc_cysc", 1);
		fChain->SetBranchStatus("dc_czsc", 1);
		fChain->SetBranchStatus("cc_segm", 1);
		fChain->SetBranchStatus("nphe", 1);
		fChain->SetBranchStatus("etot", 1);
		fChain->SetBranchStatus("ec_eo", 1);
		fChain->SetBranchStatus("ec_ei", 1);
		//added for proc_pid
		fChain->AddBranchToCache("sc_r");
		fChain->AddBranchToCache("sc_t"); 
		if (isSim){
			//mc branches
			fChain->SetBranchStatus("mcnentr",1);
			fChain->SetBranchStatus("mcid",1);
			fChain->SetBranchStatus("mcp",1);
			fChain->SetBranchStatus("mctheta",1);
			fChain->SetBranchStatus("mcphi",1);
			fChain->SetBranchStatus("mcm",1);
		}
	}
    /* *** */
    	
	treenumber = -1;

	if (isExp) {
		Info("Init", "binding e16-specific branches");
		//fChain->SetBranchAddress("evthel", &data.h10.evthel, &b_evthel); //e16
		//these branches are named the same but have different data types than user_ana-gsim version
		fChain->SetBranchAddress("id", data.h10.id, &b_id);
		fChain->SetBranchAddress("stat", data.h10.stat, &b_stat);
		fChain->SetBranchAddress("dc", data.h10.dc, &b_dc);
		fChain->SetBranchAddress("cc", data.h10.cc, &b_cc);
		fChain->SetBranchAddress("sc", data.h10.sc, &b_sc);
		fChain->SetBranchAddress("ec", data.h10.ec, &b_ec);
		fChain->SetBranchAddress("q", data.h10.q, &b_q);
		fChain->SetBranchAddress("dc_sect", data.h10.dc_sect, &b_dc_sect);
		fChain->SetBranchAddress("dc_stat", data.h10.dc_stat, &b_dc_stat);
		fChain->SetBranchAddress("ec_stat", data.h10.ec_stat, &b_ec_stat);
		fChain->SetBranchAddress("ec_sect", data.h10.ec_sect, &b_ec_sect);
		fChain->SetBranchAddress("sc_sect", data.h10.sc_sect, &b_sc_sect);
		fChain->SetBranchAddress("sc_pd", data.h10.sc_pd, &b_sc_pd);
		fChain->SetBranchAddress("sc_stat", data.h10.sc_stat, &b_sc_stat);
		fChain->SetBranchAddress("cc_sect", data.h10.cc_sect, &b_cc_sect);
		fChain->SetBranchAddress("nphe", data.h10.nphe, &b_nphe);

	}
	fChain->SetBranchAddress("evntid", &data.h10.evntid, &b_evntid);
	fChain->SetBranchAddress("npart", &data.h10.npart, &b_npart);
	fChain->SetBranchAddress("q_l", &data.h10.q_l, &b_q_l);
	fChain->SetBranchAddress("t_l", &data.h10.t_l, &b_t_l);
	fChain->SetBranchAddress("tr_time", &data.h10.tr_time, &b_tr_time);
	fChain->SetBranchAddress("gpart", &data.h10.gpart, &b_gpart);
	fChain->SetBranchAddress("p", data.h10.p, &b_p);
	fChain->SetBranchAddress("m", data.h10.m, &b_m);
	fChain->SetBranchAddress("b", data.h10.b, &b_b);
	fChain->SetBranchAddress("cx", data.h10.cx, &b_cx);
	fChain->SetBranchAddress("cy", data.h10.cy, &b_cy);
	fChain->SetBranchAddress("cz", data.h10.cz, &b_cz);
	fChain->SetBranchAddress("vx", data.h10.vx, &b_vx);
	fChain->SetBranchAddress("vy", data.h10.vy, &b_vy);
	fChain->SetBranchAddress("vz", data.h10.vz, &b_vz);
	fChain->SetBranchAddress("dc_part", &data.h10.dc_part, &b_dc_part);
	fChain->SetBranchAddress("dc_xsc", data.h10.dc_xsc, &b_dc_xsc);
	fChain->SetBranchAddress("dc_ysc", data.h10.dc_ysc, &b_dc_ysc);
	fChain->SetBranchAddress("dc_zsc", data.h10.dc_zsc, &b_dc_zsc);
	fChain->SetBranchAddress("dc_cxsc", data.h10.dc_cxsc, &b_dc_cxsc);
	fChain->SetBranchAddress("dc_cysc", data.h10.dc_cysc, &b_dc_cysc);
	fChain->SetBranchAddress("dc_czsc", data.h10.dc_czsc, &b_dc_czsc);
	fChain->SetBranchAddress("ec_part", &data.h10.ec_part, &b_ec_part);
	fChain->SetBranchAddress("etot", data.h10.etot, &b_etot);
	fChain->SetBranchAddress("ec_ei", data.h10.ec_ei, &b_ec_ei);
	fChain->SetBranchAddress("ec_eo", data.h10.ec_eo, &b_ec_eo);
	fChain->SetBranchAddress("sc_part", &data.h10.sc_part, &b_sc_part);
	fChain->SetBranchAddress("sc_t", data.h10.sc_t, &b_sc_t);
	fChain->SetBranchAddress("sc_r", data.h10.sc_r, &b_sc_r);
	fChain->SetBranchAddress("cc_part", &data.h10.cc_part, &b_cc_part);
	fChain->SetBranchAddress("cc_segm", data.h10.cc_segm, &b_cc_segm);

        /*atrivedi: following change made so that a root file with just 'thrown'
                    events sets isSim = kTRUE. This is needed so that 'anah5D' can
                    test thrown statistics */
	//if (fChain->GetBranch("l2bit")) {
        if (fChain->GetBranch("mcnentr")) {
		isSim = kTRUE;
		isExp = kFALSE;
		data.h10.is_sim = isSim;
		/* ***** ANA_USER FROM GSIM BRANCHES ***** */
		Info("Init", "binding ana_user-from-gsim-specific branches\n");
		fChain->SetBranchAddress("vidmvrt", &data.h10.vidmvrt, &b_vidmvrt);
		fChain->SetBranchAddress("ntrmvrt", &data.h10.ntrmvrt, &b_ntrmvrt);
		fChain->SetBranchAddress("xmvrt", &data.h10.xmvrt, &b_xmvrt);
		fChain->SetBranchAddress("ymvrt", &data.h10.ymvrt, &b_ymvrt);
		fChain->SetBranchAddress("zmvrt", &data.h10.zmvrt, &b_zmvrt);
		fChain->SetBranchAddress("ch2mvrt", &data.h10.ch2mvrt, &b_ch2mvrt);
		fChain->SetBranchAddress("cxxmvrt", &data.h10.cxxmvrt, &b_cxxmvrt);
		fChain->SetBranchAddress("cxymvrt", &data.h10.cxymvrt, &b_cxymvrt);
		fChain->SetBranchAddress("cxzmvrt", &data.h10.cxzmvrt, &b_cxzmvrt);
		fChain->SetBranchAddress("cyymvrt", &data.h10.cyymvrt, &b_cyymvrt);
		fChain->SetBranchAddress("cyzmvrt", &data.h10.cyzmvrt, &b_cyzmvrt);
		fChain->SetBranchAddress("stamvrt", &data.h10.stamvrt, &b_stamvrt);
		fChain->SetBranchAddress("mcnentr", &data.h10.mcnentr, &b_mcnentr);
		fChain->SetBranchAddress("mcnpart", &data.h10.mcnpart, &b_mcnpart);
		fChain->SetBranchAddress("mcst", data.h10.mcst, &b_mcst);
		fChain->SetBranchAddress("mcid", data.h10.mcid, &b_mcid);
		fChain->SetBranchAddress("mcpid", data.h10.mcpid, &b_mcpid);
		fChain->SetBranchAddress("mctheta", data.h10.mctheta, &b_mctheta);
		fChain->SetBranchAddress("mcphi", data.h10.mcphi, &b_mcphi);
		fChain->SetBranchAddress("mcp", data.h10.mcp, &b_mcp);
		fChain->SetBranchAddress("mcm", data.h10.mcm, &b_mcm);
		fChain->SetBranchAddress("mcvx", data.h10.mcvx, &b_mcvx);
		fChain->SetBranchAddress("mcvy", data.h10.mcvy, &b_mcvy);
		fChain->SetBranchAddress("mcvz", data.h10.mcvz, &b_mcvz);
		fChain->SetBranchAddress("mctof", data.h10.mctof, &b_mctof);
		//these branches are named the same but have different data types than e16 h10 tree
		fChain->SetBranchAddress("id", id, &b_id);
		fChain->SetBranchAddress("stat", stat, &b_stat);
		fChain->SetBranchAddress("dc", dc, &b_dc);
		fChain->SetBranchAddress("cc", cc, &b_cc);
		fChain->SetBranchAddress("sc", sc, &b_sc);
		fChain->SetBranchAddress("ec", ec, &b_ec);
		fChain->SetBranchAddress("q", q, &b_q);
		fChain->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
		fChain->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
		fChain->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
		fChain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
		fChain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
		fChain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
		fChain->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
		fChain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
		fChain->SetBranchAddress("nphe", nphe, &b_nphe);
	}
	
	Info("Init", "Exiting\n");
	return;
	
}

Bool_t SelH10::Notify() {
	Info("Notify", "");
		
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.
	Info("Notify", "Exiting\n");
	return kTRUE;
}

Int_t SelH10::GetEntry(Long64_t entry, Int_t getall) {
	Int_t retval = 0;
	if (fChain) {
		retval = fChain->GetTree()->GetEntry(entry, getall);
		//cout <<"Numer of bytes read for entry " << entry << " = " << retval << endl;
		Reconcile();
		//return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
	}
}

Bool_t SelH10::Process(Long64_t entry)
{
	//Info("Process", "");
	data.Clear();
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
		Int_t prevRun = data.h10.run;
		TChain *c = (TChain*) fChain;
		treenumber = c->GetTreeNumber();
		TString fn = c->GetCurrentFile()->GetName(); // was earlier GetFile() which crashed when using Prooflite
		memcpy(data.h10.fn,fn.Data(),256);
		Int_t filenumber = c->GetTreeNumber()+1; //filenumber++;
		TRegexp _regexp_run("[0-9][0-9][0-9][0-9][0-9]");
		TString rnstr = fn(_regexp_run);
		data.h10.run = rnstr.Atoi();
		Info("Process", " opened file #%d: %s (rn #%d)\n",filenumber,fn.Data(),data.h10.run);
		fFileOut->cd();
		if (strcmp(lastProcName->GetName(), "tops") != 0) fChain->GetTree()->CopyAddresses(tH10skim); //atrivedi 012013: no need to make skim after tops
		firstFileOfRun = (prevRun != data.h10.run);
	}
	//atrivedi: 09-10-12: Cofirm this logic with Evan
	//if (isExp && firstFileOfRun && data.h10.q_l == 0) return false;
	
	//atrivedi
	//tH10skim->Fill();
		
	vector<EpProcessor*>::iterator itrProcs;
	for (itrProcs = processors.begin(); itrProcs < processors.end(); itrProcs++){
		(*itrProcs)->pass = kFALSE;
	}
	top.handle(&data);
	if (processors.back()->pass){
		tOmega->Fill();
		if (strcmp(lastProcName->GetName(), "tops") != 0) tH10skim->Fill(); //atrivedi 012013: no need to make skim after tops
		nInFcBlock++;
		numfilled++;
	}
	fLumBlock.nEvnt++;
	if (data.h10idxE==0) {
		fLumBlock.nEvnt_e++;
		if (data.h10idxP>0) {
			fLumBlock.nEvnt_ep++;
			if (data.h10idxPip>0) fLumBlock.nEvnt_eppip++;
			if (data.h10idxPim>0) fLumBlock.nEvnt_eppim++;
			if (data.h10idxPip>0 && data.h10idxPim>0) fLumBlock.nEvnt_eppippim++;
		}
	}
	if (isExp && fcq_tmp != data.h10.q_l && data.h10.q_l>0) {
		fLumBlock.q_lb = data.h10.q_l-fcq_tmp;
		fLumBlock.t_l = data.h10.t_l;
		FillFc();
		fcq_tmp = data.h10.q_l;
		nfc++;
	}
	//Info("Process", "Exiting\n");
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
		//Write Objects
		if(tH10skim){
			tOmega->Write();
			tQnorm->Write();
			tLumBlock->Write();
			tH10skim->Write();
			Info("SlaveTerminate","Going to write hists from processors");
			Info("SlaveTerminate","last processor = %s", lastProcName->GetName());
			vector<EpProcessor*>::iterator itrProcs;
			for (itrProcs = processors.begin(); itrProcs < processors.end(); itrProcs++){
				Info("SlaveTerminate","Going to call ep_processor::write()");
				(*itrProcs)->write();
			}
			cleanup = kFALSE;
		}
		gDirectory = savedir;
		fFileOut->Close();
		// Cleanup or register
		if (cleanup) {
			Info("SlaveTerminate", "nothing to save: just cleanup everything ...");
			TUrl uf(*(fFileOut->GetEndpointUrl()));
			SafeDelete(fFileOut);
			gSystem->Unlink(uf.GetFile());
			SafeDelete(fProofFile);
		} else {
			Info("SlaveTerminate", "objects saved into '%s%s': sending related TProofOutputFile ...",
					fProofFile->GetFileName(), fProofFile->GetOptionsAnchor());
			fProofFile->Print();
			fOutput->Add(fProofFile);
		}
	}else {
		Info("SlaveTerminate", "Going to write non-Proof file %s", fFileOut->GetName());
		fFileOut->Write();
		fFileOut->Close();
	}
	
	Float_t ttime = swmain.RealTime();
	Float_t percentProcessed = (Float_t)eventnum/numInChain*100;
	Info("SlaveTerminate","Total: (%.2f\%) %lld/%.2f = %i events/sec",percentProcessed,eventnum,ttime,(Int_t)(eventnum/ttime));
	//while (!processors.empty()) processors.pop_back();
	//delete lastProcName;
	//delete tOmega;
	//delete tH10skim;
	//delete tQnorm;
	//delete tLumBlock;
	
	Info("SlaveTerminate", "Exiting\n");
	return;

}

void SelH10::Terminate()
{
	Info("Terminate", "");
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.
	
	//atrivedi: Following block of code moved into SlaveBegin() to make Sel_h10 work with PROOF
	/*Float_t ttime = swmain.RealTime();
	Float_t percentProcessed = (Float_t)eventnum/numInChain*100;
	printf("Total: (%.2f\%) %lld/%.2f = %i events/sec\n",percentProcessed,eventnum,ttime,(Int_t)(eventnum/ttime));
	//printf("counter = %d\n", counter);
	fFileOut->cd();
	tOmega->Write();
	tH10skim->Write();
	tQnorm->Write();
	tLumBlock->Write();
	fFileOut->Write();
	while (!processors.empty()) processors.pop_back();
	delete lastProcName;
	delete tOmega;
	delete tH10skim;
	delete tQnorm;
	delete tLumBlock;*/
	
	Info("Terminate", "Exiting\n");
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

void SelH10::FillFc()
{
	fLumBlock.run = data.h10.run;
	Int_t lastEntry = numfilled-1;
	Int_t firstEntry = lastEntry-nInFcBlock+1;
	if (firstEntry >= 0 && lastEntry > 0) {
		fLumBlock.q_l = data.h10.q_l;
		fLumBlock.evnt_first = firstEntry;
		fLumBlock.evnt_last = lastEntry;
		fLumBlock.id_last = data.h10.evntid;
		tLumBlock->Fill();
		//printf("%d - %d\t%d\n",firstEntry,lastEntry,nInFcBlock);
		for (int entry = firstEntry; entry <= lastEntry; entry++) {
			tQnorm->Fill();
		}
	}
	nInFcBlock = 0;
	fLumBlock.nEvnt = fLumBlock.nEvnt_e = fLumBlock.nEvnt_ep = fLumBlock.nEvnt_eppip
		= fLumBlock.nEvnt_eppim = fLumBlock.nEvnt_eppippim = 0;
}

void SelH10::Reconcile() {
	//Info("Reconcile", "\n");
	if (isSim) {
		//printf("Reconciling for simulated data\n");
		for (int i = 0; i < data.h10.gpart; i++) {
			data.h10.id[i] = id[i];
			//data.h10.stat[i] = stat[i];
			data.h10.dc[i] = dc[i];
			data.h10.cc[i] = cc[i];
			data.h10.sc[i] = sc[i];
			data.h10.ec[i] = ec[i];
			data.h10.q[i] = q[i];
		}
		/*for (int i = 0; i < data.h10.dc_part; i++) {
			data.h10.dc_sect[i] = dc_sect[i];
			data.h10.dc_stat[i] = dc_stat[i];
		}*/
		/*for (int i = 0; i < data.h10.ec_part; i++) {
			data.h10.ec_stat[i] = ec_stat[i];
			data.h10.ec_sect[i] = ec_sect[i];
		}*/
		for (int i = 0; i < data.h10.sc_part; i++) {
			//data.h10.sc_stat[i] = sc_stat[i];
			data.h10.sc_sect[i] = sc_sect[i];
			data.h10.sc_pd[i] = sc_pd[i];
		}
		/*for (int i = 0; i < data.h10.cc_part; i++) {
			data.h10.cc_sect[i] = cc_sect[i];
			data.h10.nphe[i] = nphe[i];
		}*/
	}
}



