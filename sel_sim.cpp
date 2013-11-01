#define SelSim_cxx
#include "sel_sim.h"
#include <TH2.h>
#include <TStyle.h>

void SelSim::Begin(TTree * /*tree*/) {

	TString option = GetOption();

}

void SelSim::SlaveBegin(TTree * /*tree*/) {
	TString option = GetOption();

}

Bool_t SelSim::Process(Long64_t entry) {

	return kTRUE;
}

void SelSim::SlaveTerminate() {
}

void SelSim::Terminate() {
}

void SelSim::Reconcile() {
	for (int i = 0; i < data.h10.gpart; i++) {
		data.h10.id[i] = id[i];
		data.h10.stat[i] = stat[i];
		data.h10.dc[i] = dc[i];
		data.h10.cc[i] = cc[i];
		data.h10.sc[i] = sc[i];
		data.h10.ec[i] = ec[i];
		data.h10.lec[i] = lec[i];
		data.h10.q[i] = q[i];
	}
	for (int i = 0; i < data.h10.dc_part; i++) {
		data.h10.dc_sect[i] = dc_sect[i];
		data.h10.dc_trk[i] = dc_trk[i];
		data.h10.dc_stat[i] = dc_stat[i];
	}
	for (int i = 0; i < data.h10.ec_part; i++) {
		data.h10.ec_stat[i] = ec_stat[i];
		data.h10.ec_sect[i] = ec_sect[i];
	}
	for (int i = 0; i < data.h10.sc_part; i++) {
		data.h10.sc_stat[i] = sc_stat[i];
		data.h10.sc_sect[i] = sc_sect[i];
		data.h10.sc_hit[i] = sc_hit[i];
		data.h10.sc_pd[i] = sc_pd[i];
	}
	for (int i = 0; i < data.h10.cc_part; i++) {
		data.h10.cc_hit[i] = cc_hit[i];
		data.h10.cc_sect[i] = cc_sect[i];
		data.h10.nphe[i] = nphe[i];
	}
}
