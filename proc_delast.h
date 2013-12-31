#ifndef PROCDELAST_H
#define PROCDELAST_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

using namespace TMath;
using namespace ParticleConstants;

class ProcDelast : public EpProcessor {

public:
	ProcDelast(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcDelast();
	void handle();
	void write();
protected:
	/*static const Int_t NUM_EVTCUTS = 16;
	enum { EVT_NULL,  EVT,       EVT_GPART0,  EVT_GPARTEQ1, EVT_GPARTEQ2, EVT_GPARTEQ3, EVT_GPARTEQ4, EVT_GPART4,
	       EVT_GOODE, EVT_GOODP, EVT_GOODPIP, EVT_GOODPIM,
	       EVT_T1,    EVT_T2,    EVT_T3,      EVT_T4,       EVT_OTHER
	     };*/
};

ProcDelast::ProcDelast(TDirectory *td,DataH10* dataH10,DataAna* dataAna) : EpProcessor(td, dataH10, dataAna) {
	dirout->cd();
	//hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	
}

ProcDelast::~ProcDelast() {
	//delete hevtsum;
	
}

void ProcDelast::handle() {
	//Info("In ProcDelast::handle()");
	pass = kFALSE;
	
	if (pass) {
		EpProcessor::handle();
	}
}




void ProcDelast::write(){

}

#endif // PROCDELAST_H
