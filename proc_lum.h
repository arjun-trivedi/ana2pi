/******************************
/* The following can be used as a template for
/* processors that make data, for example:
/* proc_delast,
/* proc_d2pi 
/* etc 
********************************/

#ifndef PROCLUM_H
#define PROCLUM_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include "proc_eid.h"

using namespace TMath;
using namespace ParticleConstants;

class ProcLum : public EpProcessor {

public:
	ProcLum(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcLum();
	void handle();
	void write();
protected:
	ProcEid* proc_eid;
	static const Int_t NUM_EVTCUTS = 2;
	enum {EVT_NULL,EVT,EVT_GOODE};
	TH1F* _hQ;
	float _totalQ;
 	float _qcurr;
 	float _qprev;
 	float _deltaq;
 	float _q_temp;
};

ProcLum::ProcLum(TDirectory *td,DataH10* dataH10,DataAna* dataAna) : EpProcessor(td, dataH10, dataAna) {
	proc_eid=new ProcEid(td->mkdir("eid"),dataH10,dataAna);
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	_hQ=new TH1F("hQ","Total charge",1,0,1);
	_totalQ=0.0;
 	_qcurr=0.0 ;
 	_qprev=0.0 ;
 	_deltaq=0.0 ;
 	_q_temp=0.0;
	
}

ProcLum::~ProcLum() {
	//delete hevtsum;
	
}

void ProcLum::handle() {
	//Info("In ProcLum::handle()","");
	hevtsum->Fill(EVT);

	_q_temp=dH10->q_l;
	_qcurr=_q_temp;
	//cout<<"q_l="<<q_l<<endl;

	if ((_q_temp>0.)){
		// cout<<"q_l"<<q_temp<<"qcurr"<<qcurr<<endl;
		if(( _qcurr>_qprev && _qprev!= 0)){
			_deltaq = _qcurr-_qprev;
			_totalQ += _deltaq;
			//cout<<"qcurr="<<qcurr<<"qprev="<<qprev<<"deltaq="<<deltaq<<"totalQ="<<totalQ<<endl;
			//printf("qcurr=%.2f,qprev=%.2f,deltaq=%.2f,totalQ=%.2f",_qcurr,_qprev,_deltaq,_totalQ);
		// cout<<"q_l"<<q_temp<<"qcurr"<<qcurr<<"qprev"<<qprev<<"deltaq"<<deltaq<<endl;
		}
	_qprev = _qcurr;
	_hQ->SetBinContent(1,_totalQ);
	//cout<<"totalQ="<<totalQ<<endl;
    }

	pass = kFALSE;
	Bool_t gE = kFALSE;
	gE=proc_eid->goodE();
	if (gE) {
		hevtsum->Fill(EVT_GOODE);
		EpProcessor::handle();
	}
}




void ProcLum::write(){

}

#endif // PROCLUM_H
