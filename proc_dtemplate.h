/******************************
/* The following can be used as a template for
/* processors that make data, for example:
/* proc_delast,
/* proc_d2pi 
/* etc 
********************************/

#ifndef PROCD<TEMPLATE>_H
#define PROCD<TEMPLATE>_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

using namespace TMath;
using namespace ParticleConstants;

class ProcD<template> : public EpProcessor {

public:
	ProcD<template>(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcD<template>();
	void handle();
	void write();
protected:
	/*static const Int_t NUM_EVTCUTS = 16;
	enum { EVT_NULL,  EVT,       EVT_GPART0,  EVT_GPARTEQ1, EVT_GPARTEQ2, EVT_GPARTEQ3, EVT_GPARTEQ4, EVT_GPART4,
	       EVT_GOODE, EVT_GOODP, EVT_GOODPIP, EVT_GOODPIM,
	       EVT_T1,    EVT_T2,    EVT_T3,      EVT_T4,       EVT_OTHER
	     };*/
};

ProcD<template>::ProcD<template>(TDirectory *td,DataH10* dataH10,DataAna* dataAna) : EpProcessor(td, dataH10, dataAna) {
	dirout->cd();
	//hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	
}

ProcD<template>::~ProcD<template>() {
	//delete hevtsum;
	
}

void ProcD<template>::handle() {
	//Info("In ProcD<template>::handle()");
	pass = kFALSE;
	
	if (pass) {
		EpProcessor::handle();
	}
}




void ProcD<template>::write(){

}

#endif // PROCD<TEMPLATE>_H
