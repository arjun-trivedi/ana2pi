//at-h8 #include "q2w_bng.h"
#include "data_h10.h"
#include "data_ana.h"
#include "h10looper.h"
#include "ep_processor.h"
#include "proc_eid.h"
#include "proc_delast.h"
#include "proc_efid.h"
#include "proc_pid.h"
#include "proc_skim_q.h"
//at-h8 #include "proc_skim_q2w.h"
//at-h8 #include "proc_elist_q2w.h"
#include "proc_mom_cor.h"
#include "proc_d2pi.h"
#include "proc_copy_h10.h"

#include <TChain.h>
#include <TFileCollection.h>
#include <TFile.h>
#include <TRegexp.h>

#include <unistd.h>
#include <stdio.h>

using namespace std;

//user input
TString fin=""; //h10.lst
TString procorder="";
TString output="";
/*TString expt="";
TString dtyp="";
TString rctn="";*/ 
TString h10type="";

TString str_use_q2w_elist;
Bool_t use_q2w_elist=kFALSE;
TString q2w="";
TString fname_q2w_el="";
TFile* f_q2w_el=NULL;

TString str_nentries="";
Long64_t nentries=1000000000;

//objects setup by proc_h10
//TString h10type;
TChain* h10chain;
TFile* fout;
DataH10* dH10;
DataAna* dAna;
EpProcessor* proc_chain;
H10Looper* h10looper;

void parseArgs(int argc, char* const argv[]);
TDirectory* mkdir(const char* dirname);
EpProcessor* SetupProcs();
EpProcessor* SetupProcSkimQ2W(TString proc_q2wskim_dcptr);

int main(int argc,  char* const argv[])
{
	Info("proc_h10::main()", "\n");
	parseArgs(argc, argv);
	
	/*if (fin==""||procorder==""||expt==""||dtyp==""||rctn==""){*/
	if (fin==""||h10type==""||procorder==""||output==""){
		printf("Not all arguments entered\n");
		printf("Execute \"proc_h10 -h\" to see correct usage\n");
		return 0;
	}

	TObjArray *h10type_tokens = h10type.Tokenize(":");
	TString expt  = h10type_tokens->At(0)->GetName();
	TString dtyp = h10type_tokens->At(1)->GetName();
	TString rctn = h10type_tokens->At(2)->GetName();

	if (expt!="e1f" && expt!="e16"){
		printf("Incorrect expt entered: %s\n", expt.Data());
		return 0;
	}
	if (dtyp!="exp" && dtyp!="sim"){
		printf("Incorrect dtyp entered: %s\n", dtyp.Data());
		return 0;
	}
	if (rctn!="2pi" && rctn!="elast" && rctn!="2pi_userana" && rctn!="elast_userana"){
		printf("Incorrect rctn entered: %s\n", rctn.Data());
		return 0;
	}
	if (use_q2w_elist){
		Info("proc_h10","Going to use q2w_elist\n");
		f_q2w_el=new TFile(fname_q2w_el);
		// printf("%s",f_q2w_el->GetName());
	}else{
		Info("proc_h10","Not going to use q2w_elist\n");	
	}


	h10chain = new TChain("h10");
	TFileCollection fc("fileList", "", fin.Data());
	h10chain->AddFileInfoList((TCollection*) fc.GetList());
		
	//TString fout_name = TString::Format("d%s.root",rctn.Data());
	//TString fout_name = TString::Format("d%s.root",rctn.Data());
	
	fout = new TFile(output,"RECREATE");

	h10type=TString::Format("%s:%s:%s",expt.Data(),dtyp.Data(),rctn.Data());
	dH10 = new DataH10(h10type);
	dAna = new DataAna();
	proc_chain = SetupProcs();
	h10looper = new H10Looper(h10chain,dH10,dAna,proc_chain,use_q2w_elist,f_q2w_el,q2w);
	h10looper->Loop(nentries);

	fout->Write();
	fout->Close();
	
	return 0;
}

void parseArgs(int argc, char* const argv[]){
	int c;
	//char *filename;
	extern char *optarg;
	extern int optind, optopt, opterr;

	while ((c = getopt(argc, argv, "hi:t:p:o:l:m:q:n:")) != -1) {
		switch(c) {
		case 'h':
			printf("proc_h10 -i <h10.lst> -t <expt>:<dtyp>:<rctn> -p <procorder> -o <output> -l <use_q2w_elist> -m <fname_q2w_el> -q=<q2wX>-n <nevts>\n");
			printf("<expt>=e1f/e16\n");
			printf("<dtyp>=exp/sim\n");
			printf("<rctn>=2pi/elast/2pi_userana/elast_userana\n");
			printf("<use_q2w_elist>=True/[false]\n");
			break;
		case 'i':
			fin = optarg;
			break;	
		case 't':
			h10type = optarg;
			break;
		case 'p':
			procorder = optarg;
			break;
		case 'o':
			output = optarg;
			break;
		case 'l':
			str_use_q2w_elist = optarg;
			if (str_use_q2w_elist.EqualTo("True")){
				use_q2w_elist=kTRUE;
			}
			break;
		case 'm':
			fname_q2w_el=optarg;
			break;
		case 'q':
			q2w=optarg;
			break;
		case 'n':
			str_nentries = optarg;
			nentries = str_nentries.Atoll();
			break;
		case ':':
			printf("-%c without options\n", optopt);
			break;
		case '?':
			printf("unknown arg %c\n", optopt);
			break;
		}
	}
	return;
}


TDirectory* mkdir(const char* dirname)
{
	TString dname = dirname;
	TDirectory *ret = 0;
	Int_t numTries = 0;
	while ( ret == 0 && numTries++<10 ) {
		if (numTries>1) dname += numTries;
		ret = fout->mkdir(dname.Data());
	}
	return ret;
}

EpProcessor* SetupProcs(){
   //Get list of processors
   TCollection *procorder_tokens = procorder.Tokenize(":");
   Info("proc_h10::SetupProcs()", "procorder = %s\n", procorder.Data());
   procorder_tokens->Print();
         
   //instantiate topProc; by design, topProc "builds" a cascaded chain of Procs
   EpProcessor* proc_chain = new EpProcessor();
            
   if (procorder_tokens->IsEmpty()) {
      Info("SetupProcs", "no processors specified; building default pipeline\n");
      EpProcessor *proc;
      proc = new ProcEid(mkdir("eid"),dH10,dAna);
      //procs.push_back(proc);
      proc_chain->add(proc);
      //lastProcName= new TObjString("eid");
      //Info("SlaveBegin","last processor = %s", lastProcName->GetName());
   } else {
      TIter iter(procorder_tokens);
      while(TObjString *obj_str = (TObjString*)iter.Next()) {
         EpProcessor *proc;
         TString str = obj_str->GetString();
         if (str.EqualTo("eid"))             proc = new ProcEid(mkdir("eid"),dH10,dAna);
         else if (str.EqualTo("eidmon"))     proc = new ProcEid(mkdir("eid"),dH10,dAna,kTRUE);
         else if (str.EqualTo("eidmononly")) proc = new ProcEid(mkdir("eid"),dH10,dAna,kTRUE,kTRUE);
         else if (str.EqualTo("efid"))       proc = new ProcEFid(mkdir("fid"),dH10,dAna);
         else if (str.EqualTo("delast"))     proc = new ProcDelast(mkdir("delast"),dH10,dAna);
         else if (str.EqualTo("efidmon"))    proc = new ProcEFid(mkdir("fid"),dH10,dAna,kTRUE);
         else if (str.EqualTo("efidmononly"))proc = new ProcEFid(mkdir("fid"),dH10,dAna,kTRUE,kTRUE);
         else if (str.EqualTo("qskim"))      proc = new ProcSkimQ(mkdir("qskim"),dH10,dAna);
         //at-h8 else if (str.EqualTo("q2wskim"))    proc = new ProcSkimQ2W(mkdir("q2wskim"),dH10,dAna);
         //at-h8 else if (str.Contains(TRegexp("q2wskim[0-9]+[0-9]?"))) proc = SetupProcSkimQ2W(str);
         //at-h8 else if (str.EqualTo("q2welist"))   proc = new ProcEListQ2W(mkdir("q2welist"),dH10,dAna);
         else if (str.EqualTo("mom"))        proc = new ProcMomCor(mkdir("mom"),dH10,dAna);
         else if (str.EqualTo("pid"))        proc = new ProcPid(mkdir("pid"),dH10,dAna);
         else if (str.EqualTo("pidmon"))     proc = new ProcPid(mkdir("pid"),dH10,dAna,kTRUE);
         else if (str.EqualTo("pidmononly")) proc = new ProcPid(mkdir("pid"),dH10,dAna,kTRUE,kTRUE);
         else if (str.EqualTo("copyh10"))    proc = new ProcCopyH10(mkdir("copyh10"),dH10,dAna);
         else if (str.EqualTo("d2piT"))      proc = new ProcD2pi(mkdir("d2piT"),dH10,dAna,kTRUE,kFALSE);
         else if (str.EqualTo("d2piR"))      proc = new ProcD2pi(mkdir("d2piR"),dH10,dAna,kFALSE,kTRUE);
         else if (str.EqualTo("d2piTR"))     proc = new ProcD2pi(mkdir("d2piTR"),dH10,dAna,kTRUE,kTRUE);
         else if (str.EqualTo("d2piT_tree"))      proc = new ProcD2pi(mkdir("d2piT"),dH10,dAna,kTRUE,kFALSE,kTRUE);
         else if (str.EqualTo("d2piR_tree"))      proc = new ProcD2pi(mkdir("d2piR"),dH10,dAna,kFALSE,kTRUE,kTRUE);
         else if (str.EqualTo("d2piTR_tree"))     proc = new ProcD2pi(mkdir("d2piTR"),dH10,dAna,kTRUE,kTRUE,kTRUE);
         else {
            Info("Init","%s unrecognized processor\n",str.Data());
            continue;
         }
         //procs.push_back(proc);
         proc_chain->add(proc);
         //lastProcName = obj_str; //works because after last iteration, lastProcName is not updated
      }
      //Info("SlaveBegin","last processor = %s", lastProcName->GetName());
   }
   return proc_chain;
}

//at-h8
/*EpProcessor* SetupProcSkimQ2W(TString proc_q2wskim_dcptr){
	Float_t q2min,q2max,wmin,wmax=0.0;
	TString str_q2w_binnum=proc_q2wskim_dcptr.Remove(0,7);//removing "q2wskim" from "q2wskimXX"
	Int_t binnum=str_q2w_binnum.Atoi();
	q2min=kQ2W_CrsBin[binnum-1].q2min;
    q2max=kQ2W_CrsBin[binnum-1].q2max;
    wmin=kQ2W_CrsBin[binnum-1].wmin;
    wmax=kQ2W_CrsBin[binnum-1].wmax;
    EpProcessor* proc = new ProcSkimQ2W(mkdir("q2wskim"),dH10,dAna,q2min,q2max,wmin,wmax);
    return proc;
}*/




