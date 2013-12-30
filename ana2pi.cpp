#include "data_h10.h"
#include "data_ana.h"
#include "ep_processor.h"
#include "proc_eid.h"
/*#include "proc_efid.h"
#include "proc_pid.h"
#include "proc_skim_q.h"
#include "proc_mom_cor.h"
#include "proc_top.h"
#include "proc_skim_q2w.h"
#include "proc_fill_skim.h"
#include "proc_copy_h10.h"*/

#include <TChain.h>
#include <TFileCollection.h>
#include <TFile.h>

#include <unistd.h>
#include <stdio.h>

using namespace std;

bool make_delast;

//user input
TString fin; //h10 or h10.lst
TString fout_name; //name of output file
TString h10type;
TString procorder;

//objects setup by ana2pi
TFile* fout;
DataH10* dH10;
DataAna* dAna;
EpProcessor* top_proc;

void parseArgs(int argc, char* const argv[]);
EpProcessor* SetupProcs();



int main(int argc,  char* const argv[])
{
	Info("ana2pi::main()", "\n");
	parseArgs(argc, argv);
	
	TChain *chain = new TChain("h10");
	if (make_delast){
		h10type="e1f:sim:elas";
		//TFileCollection fc("fileList", "", inflst.Data());
		//chain->AddFileInfoList((TCollection*) fc.GetList());
	} 
	fout = new TFile(fout_name,"RECREATE");
	dH10 = new DataH10(h10type);
	dAna = new DataAna();
	top_proc = SetupProcs();
			
	//delete chain;
	
	return 0;
}

void parseArgs(int argc, char* const argv[]){
	int c;
	//char *filename;
	extern char *optarg;
	extern int optind, optopt, opterr;

	while ((c = getopt(argc, argv, ":ei:o:p:")) != -1) {
		switch(c) {
		case 'e':
			printf("Going to make DElast\n");
			make_delast = true;
			break;
		case 'p':
			procorder = optarg;
			printf("Processor list = %s", procorder.Data());
			break;
		case 'i':
			fin = optarg;
			printf("input filename is %s\n", fin.Data());
			break;
		case 'o':
			fout_name = optarg;
			printf("output filename is %s\n", fout_name.Data());
			break;
		case ':':
			printf("-%c without filename\n", optopt);
			break;
		case '?':
			printf("unknown arg %c\n", optopt);
			break;
		}
	}
	return;
}

EpProcessor* SetupProcs(){
   //Get list of processors
   TCollection *procorder_tokens = procorder.Tokenize(":");
   Info("ana2pi::SetupProcs()", "List of processors received:\n");
   procorder_tokens->Print();
         
   //instantiate topProc; by design, topProc "builds" a cascaded chain of Procs
   EpProcessor* top_proc = new EpProcessor();
            
   if (procorder_tokens->IsEmpty()) {
      Info("SetupProcs", "no processors specified; building default pipeline\n");
      EpProcessor *proc;
      proc = new ProcEid(fout->mkdir("eid"),dH10,dAna);
      //procs.push_back(proc);
      top_proc->add(proc);
      //lastProcName= new TObjString("eid");
      //Info("SlaveBegin","last processor = %s", lastProcName->GetName());
   } else {
      TIter iter(procorder_tokens);
      while(TObjString *obj_str = (TObjString*)iter.Next()) {
         EpProcessor *proc;
         TString str = obj_str->GetString();
         if (str.EqualTo("eid"))            proc = new ProcEid(fout->mkdir("eid"),dH10,dAna);
         else if (str.EqualTo("eidmon"))     proc = new ProcEid(fout->mkdir("eid"),dH10,dAna,kTRUE);
         else if (str.EqualTo("eidmononly")) proc = new ProcEid(fout->mkdir("eid"),dH10,dAna,kTRUE,kTRUE);
         /*else if (str.EqualTo("efid"))       proc = new ProcEFid(mkdir("fid"),dH10,dAna);
         else if (str.EqualTo("efidmon"))    proc = new ProcEFid(mkdir("fid"),dH10,dAna,kTRUE);
         else if (str.EqualTo("efidmononly"))proc = new ProcEFid(mkdir("fid"),dH10,dAna,kTRUE,kTRUE);
         else if (str.EqualTo("qskim"))       proc = new ProcSkimQ(mkdir("qskim"),dH10,dAna);
         else if (str.EqualTo("mom"))      proc = new ProcMomCor(mkdir("mom"),dH10,dAna);
         else if (str.EqualTo("pid"))      proc = new ProcPid(mkdir("pid"),dH10,dAna);
         else if (str.EqualTo("pidmon"))     proc = new ProcPid(mkdir("pid"),dH10,dAna,kTRUE);
         else if (str.EqualTo("pidmononly")) proc = new ProcPid(mkdir("pid"),dH10,dAna,kTRUE,kTRUE);
         else if (str.EqualTo("top"))      proc = new ProcTop(mkdir("top"),dH10,dAna);
         else if (str.EqualTo("q2wskim")) proc = new ProcSkimQ2W(mkdir("q2wskim"),dH10,dAna);
         else if (str.EqualTo("fillskim"))   proc = new ProcFillSkim(mkdir("skim"),dH10,dAna);
         else if (str.EqualTo("copyh10")) proc = new ProcCopyH10(fFileOut,dH10,dAna);*/
         else {
            Info("Init","%s unrecognized processor\n",str.Data());
            continue;
         }
         //procs.push_back(proc);
         top_proc->add(proc);
         //lastProcName = obj_str; //works because after last iteration, lastProcName is not updated
      }
      //Info("SlaveBegin","last processor = %s", lastProcName->GetName());
   }
   return top_proc;
}




