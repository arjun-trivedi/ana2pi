#include "sel_h10.h"

#include <TChain.h>
#include <TFileCollection.h>

#include <unistd.h>
#include <stdio.h>

using namespace std;

void parseArgs(int argc, char* const argv[]);

bool _useproof = false; //hard coded for now
bool _makeqskim = false;
bool _makeyields = false;

bool is_h10_Exp=false;
bool is_h10_Sim=false;
bool is_h10_SimMcOnly=false;

TString h10type;
TString inflst;
TString outfname;

int main(int argc,  char* const argv[])
{
	cout << "analysis2pi" << endl;
	parseArgs(argc, argv);
	
	TChain *chain = new TChain("h10");
	if (_makeyields) chain->Add(inflst.Data());
	else {
		TFileCollection fc("fileList", "", inflst.Data());
		chain->AddFileInfoList((TCollection*) fc.GetList());
	}
		
	SelH10 *s = new SelH10();
	if(_useproof){
		TList* inputList = new TList();
		inputList->Add(new TNamed("PROOF_OUTPUTFILE", outfname.Data()));
		//atrivedi 071313: currently setup to make qskim for e1f.exp only
		inputList->Add(new TNamed("H10TYPE", "e1f:exp"));
		s->SetInputList(inputList);
	}else{
		TList* inputList = new TList();
		inputList->Add(new TNamed("OUTPUTFILE", outfname.Data()));
		//atrivedi 071313: currently setup to make qskim for e1f.exp only
		inputList->Add(new TNamed("H10TYPE", "e1f:exp"));
		s->SetInputList(inputList);
	}
	if(_useproof){
		chain->SetProof();
	}
	
	if(_makeqskim){
		printf("Going to make qskim: %s ---> %s\n", inflst.Data(), outfname.Data());
		chain->Process(s,"eid:qskim:copyh10");
		printf("Done making qskim: %s ---> %s\n", inflst.Data(), outfname.Data());
	}
	if(_makeyields){
		if (strcmp(h10type,"exp")==0) is_h10_Exp=true;
		else if (strcmp(h10type,"sim")==0) is_h10_Sim=true;
		else if (strcmp(h10type,"simMcOnly")==0) is_h10_SimMcOnly=true;
		else {
			printf("Cannot Proceed: -y option requires -t = (exp/sim/simMcOnly)\n");
			delete chain;
			delete s;
			return 0;
		}
		if (is_h10_Exp){
			printf("Going to make yields from EXP data: %s ---> %s\n", inflst.Data(), outfname.Data());
			chain->Process(s,"eid:qskim:tops");
			printf("Done making yields from EXP data: %s ---> %s\n", inflst.Data(), outfname.Data());
		}else if(is_h10_Sim){
			printf("Going to make yields from SIM data: %s ---> %s\n", inflst.Data(), outfname.Data());
			chain->Process(s,"tops:eid:qskim:tops");
			printf("Done making yields from SIM data: %s ---> %s\n", inflst.Data(), outfname.Data());
		}else if(is_h10_SimMcOnly){
			printf("Going to make yields from MC data: %s ---> %s\n", inflst.Data(), outfname.Data());
			chain->Process(s,"tops");
			printf("Done making yields from MC data: %s ---> %s\n", inflst.Data(), outfname.Data());
		}
	}
	
	delete chain;
	delete s;
	
	return 0;
}

void parseArgs(int argc, char* const argv[]){
	int c;
	//char *filename;
	extern char *optarg;
	extern int optind, optopt, opterr;

	while ((c = getopt(argc, argv, ":qyt:i:o:")) != -1) {
		switch(c) {
		case 'q':
			printf("Going to make qskim\n");
			_makeqskim = true;
			break;
		case 'y':
			printf("Going to make yields\n");
			_makeyields = true;
			break;	
		case 't':
			h10type = optarg;
			printf("h10type is %s\n", h10type.Data());
			break;
		case 'i':
			inflst = optarg;
			printf("input filename is %s\n", inflst.Data());
			break;
		case 'o':
			outfname = optarg;
			printf("output filename is %s\n", outfname.Data());
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


