#ifndef PROCSKIMQDELAST_H
#define PROCSKIMQDELAST_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TH1.h>

using namespace ParticleConstants;

/******************************************************
[03-08-15]
+ Pass Events only if:
	1. The 1st particle is an Electron, AND
	2. The total number of Reconstructed tracks is < 4, AND
	3. The charges of the Reconstructed tracks are identified as:
		+ 1POS+0NEG+0NEUTRAL
*******************************************************/


class ProcSkimQElast : public EpProcessor
{

public:
	ProcSkimQElast(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcSkimQElast();
	
	void handle();
	
protected:
	static const Int_t NUM_EVTCUTS=4;
	enum { EVT_NULL, EVT, EVT_1ST_E, EVT_LT_EQ_4, EVT_1POS_EX};
};

ProcSkimQElast::ProcSkimQElast(TDirectory *td,DataH10* dataH10,DataAna* dataAna) 
					:EpProcessor(td, dataH10, dataAna)
{
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1ST_E,"e^{-} 1^{st}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_LT_EQ_4,"gpart <= 4");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1POS_EX,"X^{+}");
}

ProcSkimQElast::~ProcSkimQElast()
{
	delete hevtsum;
	delete hists;
}

void ProcSkimQElast::handle()
{
	//Info("ProcSkimQElast::handle()", "");
	pass = kFALSE;
		
	hevtsum->Fill(EVT);
	if ( dH10->id[0]==ELECTRON ) {
		hevtsum->Fill(EVT_1ST_E);
		if ( dH10->gpart<=4) { 
		    hevtsum->Fill(EVT_LT_EQ_4);
			Int_t numPos=0;
			Int_t numNeg=0;
			Int_t num0=0;
			for (Int_t i=1;i<dH10->gpart;i++) {
			    if      (dH10->q[i]==1)  numPos++;
			    else if (dH10->q[i]==-1) numNeg++;
			    else if (dH10->q[i]==0)  num0++;
			}
			Int_t binNum = 0;
			if (numPos==1 && numNeg==0 && num0==0) {	
				binNum=hevtsum->Fill(EVT_1POS_EX);
				dAna->skimq_elast.isEVT_1POS_EX=kTRUE;
			}			
			if (binNum>0) {
				pass=kTRUE;
				EpProcessor::handle();
			}
		} 
	}
}

#endif // PROCSKIMQDELAST_H