#include "h10looper_e1f.h"

#include <TChain.h>
#include <TFileCollection.h>
#include <TFile.h>
#include <TRegexp.h>
#include <TLorentzVector.h> //! Makes Info() work

#include <unistd.h>
#include <stdio.h>

using namespace std;

//! user input
TString h10lst="";
TString h10type="";
TString fout_name="";
//! + nentries can be specified by user
//! + default=0, in which case number of entries in TChain is used, as per logic in h10looper.C:
//!   'nentries==0?nentries_to_proc=nentries_chain:nentries_to_proc=nentries;'
TString str_nentries="";
Long64_t nentries=0;//[08-09-15]used to be 1B, till I cleaned up the "logic"

void parseArgs(int argc, char* const argv[]);

int main(int argc,  char* const argv[])
{
	Info("proc_h10_lite::main()", "\n");

	//! Get user input
	parseArgs(argc, argv);
	
	//! Make sure user input is OK
	if (h10lst==""||h10type==""||fout_name==""){
		printf("Not all arguments entered\n");
		printf("Execute \"proc_h10_lite -h\" to see correct usage\n");
		return 0;
	}
	//! Make sure h10type is OK
	TObjArray *h10type_tokens = h10type.Tokenize(":");
	TString expt  = h10type_tokens->At(0)->GetName();
	TString dtyp = h10type_tokens->At(1)->GetName();
	TString seq = h10type_tokens->At(2)->GetName();
	if (expt!="e1f" && expt!="e16"){
		printf("Incorrect expt entered: %s\n", expt.Data());
		return 0;
	}
	if (dtyp!="exp" && dtyp!="sim"){
		printf("Incorrect dtyp entered: %s\n", dtyp.Data());
		return 0;
	}
	if (seq!="recon" && seq!="thrown"){
		printf("Incorrect seq entered: %s\n", dtyp.Data());
		return 0;
	}
	Info("proc_h10_lite","h10lst=%s,h10type=%s,fout_name=%s,nentries=%llu\n",
		h10lst.Data(),h10type.Data(),fout_name.Data(),nentries);

	//! Now setup h10looper and call h10looper::Loop()
	//! h10chain
	TChain*	h10chain = new TChain("h10");
	TFileCollection fc("fileList", "", h10lst.Data());
	h10chain->AddFileInfoList((TCollection*) fc.GetList());
	//! Set up h10looper	
	h10looper_e1f* h10lpr_e1f=new h10looper_e1f(h10chain,h10type,fout_name,nentries);
	//h10looper.Loop()
	

	delete h10chain;
	delete h10lpr_e1f;
	
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
			printf("proc_h10_lite -i <h10.lst> -t <expt>:<dtyp>:<seq> -o <fout_name> -n <nevts>\n");
			printf("<expt>=e1f/e16\n");
			printf("<dtyp>=exp/sim\n");
			printf("<seq>=recon/thrown\n");
			break;
		case 'i':
			h10lst = optarg;
			break;	
		case 't':
			h10type = optarg;
			break;
		case 'o':
			fout_name = optarg;
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
