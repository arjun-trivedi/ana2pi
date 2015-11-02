#include "h10looper_e1f.h"
#include "h10looper_2pi.h"

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
TString cutsncors="";
TString fout_name="";
//! + nentries can be specified by user
//! + default=0, in which case number of entries in TChain is used, as per logic in h10looper.C:
//!   'nentries==0?nentries_to_proc=nentries_chain:nentries_to_proc=nentries;'
TString str_nentries="";
Long64_t nentries=0;//[08-09-15]used to be 1B, till I cleaned up the "logic"

//! cuts to be made in addition to 'dflt'
//! syntax= 1:2:3:4: i.e. number followed by colon
TString adtnl_opts="";

void parseArgs(int argc, char* const argv[]);

int main(int argc,  char* const argv[])
{
	Info("proc_h10_lite::main()", "\n");

	//! Get user input
	parseArgs(argc, argv);
	
	//! Make sure user input is OK
	if (h10lst==""||h10type==""||cutsncors==""||fout_name==""){
		printf("Not all arguments entered\n");
		printf("Execute \"proc_h10_lite -h\" to see correct usage\n");
		return 0;
	}
	//! Make sure h10type is OK
	TObjArray *h10type_tokens = h10type.Tokenize(":");
	TString expt  = h10type_tokens->At(0)->GetName();
	TString dtyp = h10type_tokens->At(1)->GetName();
	TString rctn = h10type_tokens->At(2)->GetName();
	TString seq = h10type_tokens->At(3)->GetName();
	if (expt!="e1f" && expt!="e16"){
		printf("Incorrect expt entered: %s\n", expt.Data());
		return 0;
	}
	if (dtyp!="exp" && dtyp!="sim"){
		printf("Incorrect dtyp entered: %s\n", dtyp.Data());
		return 0;
	}
	if (rctn!="2pi" && rctn!="elast"){
		printf("Incorrect rctn entered: %s\n", rctn.Data());
		return 0;
	}
	if (seq!="recon" && seq!="thrown"){
		printf("Incorrect seq entered: %s\n", dtyp.Data());
		return 0;
	}
	Info("proc_h10_lite","h10lst=%s,h10type=%s,cutsncors=%s,fout_name=%s,nentries=%llu\n,adtnl_opts=%s",
		h10lst.Data(),h10type.Data(),cutsncors.Data(),fout_name.Data(),nentries,adtnl_opts.Data());

	//! Now setup h10looper and call h10looper::Loop()
	//! h10chain
	TChain*	h10chain = new TChain("h10");
	TFileCollection fc("fileList", "", h10lst.Data());
	h10chain->AddFileInfoList((TCollection*) fc.GetList());
	//! Set up h10looper	
	if (rctn=="2pi"){
		h10looper_2pi* h10lpr_2pi=new h10looper_2pi(h10type,h10chain,cutsncors,fout_name,nentries,adtnl_opts);
		h10lpr_2pi->Loop();
		delete h10chain;
		delete h10lpr_2pi;
	}else if (rctn=="elast"){
		h10looper_e1f* h10lpr_e1f=new h10looper_e1f(h10type,h10chain,cutsncors,fout_name,nentries,adtnl_opts);
		h10lpr_e1f->Loop();	
		delete h10chain;
		delete h10lpr_e1f;
	}
	
	return 0;
}

/*
Following function implemented as per:
http://www.gnu.org/software/libc/manual/html_node/Example-of-Getopt.html
*/
void parseArgs(int argc, char* const argv[]){
	int c;
	extern char *optarg;
	extern int optind, optopt, opterr;
	int index;

	while ((c = getopt(argc, argv, "hi:t:c:o:n:")) != -1) {
		switch(c) {
		case 'h':
			printf("proc_h10_lite -i <h10.lst> -t <expt>:<dtyp>:<rctn>:<seq> -c <cutsncors> -o <fout_name> -n <nevts> [adtnl_opts]\n");
			printf("<expt>=e1f/e16\n");
			printf("<dtyp>=exp/sim\n");
			printf("<rcnt>=2pi/elast\n");
			printf("<seq>=recon/thrown\n");
			break;
		case 'i':
			h10lst = optarg;
			break;	
		case 't':
			h10type = optarg;
		case 'c':
			cutsncors = optarg;
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

	//! Use the following top debug argc,argv
	/*printf("argc=%d\n",argc);
	printf("argv[0]=%s\n",argv[0]);
	printf("argv[1]=%s\n",argv[1]);
	printf("argv[%d]=%s\n",argc-1,argv[argc-1]);
	printf("optind=%d\n",optind);*/

	//! Get adtnl_opts, which is the last option passed to the program
	adtnl_opts=argv[argc-1];
	/*for (index=optind;index<argc;index++){
		string arg=argv[index];
	}*/
  	
	return;
}
