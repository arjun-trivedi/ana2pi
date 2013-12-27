#include "data_h10.h"

#include <TChain.h>
#include <TFileCollection.h>

#include <unistd.h>
#include <stdio.h>

using namespace std;

void parseArgs(int argc, char* const argv[]);

bool make_delast;

TString h10type;
TString inflst;
TString outfname;

int main(int argc,  char* const argv[])
{
	Info("ana2pi::main()", "Could not determine h10typ.exp!\n");
	parseArgs(argc, argv);
	
	TChain *chain = new TChain("h10");
	if (make_delast){
		h10type="e1f:sim:elas";
		TFileCollection fc("fileList", "", inflst.Data());
		chain->AddFileInfoList((TCollection*) fc.GetList());
	} 

	DataH10* dH10 = new DataH10(h10type);
			
	delete chain;
	
	return 0;
}

void parseArgs(int argc, char* const argv[]){
	int c;
	//char *filename;
	extern char *optarg;
	extern int optind, optopt, opterr;

	while ((c = getopt(argc, argv, ":e:i:o:")) != -1) {
		switch(c) {
		case 'e':
			printf("Going to make DElast\n");
			make_delast = true;
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


