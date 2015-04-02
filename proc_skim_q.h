#ifndef PROCSKIMQ_H
#define PROCSKIMQ_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TH1.h>

using namespace ParticleConstants;

/******************************************************
[02-24-15]
+ Pass Events only if:
	1. The 1st particle is an Electron, AND
	2. The charges of the Reconstructed tracks are identified as either:
		+ At least 2POS AND At least 1NEG (top1 candidate)
			 OR
		+ At least 2POS (top2 candidate)
			 OR
		+ At least 1POS1NEG (top3 candidate or top4 candidate)
*******************************************************/


class ProcSkimQ : public EpProcessor
{

public:
	ProcSkimQ(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcSkimQ();
	
	void handle();
	
protected:
	static const Int_t NUM_EVTCUTS=5;
	enum { EVT_NULL, EVT, EVT_1ST_E, EVT_ETGT_2POS_ETGT_1NEG, EVT_ETGT_2POS, EVT_ETGT_1POS_ETGT_1NEG };
};

ProcSkimQ::ProcSkimQ(TDirectory *td,DataH10* dataH10,DataAna* dataAna) 
					:EpProcessor(td, dataH10, dataAna)
{
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1ST_E,"e^{-} 1^{st}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ETGT_2POS_ETGT_1NEG,"X^{+}>=2 + X^{-}>=1");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ETGT_2POS,"X^{+}>=2");
	hevtsum->GetXaxis()->SetBinLabel(EVT_ETGT_1POS_ETGT_1NEG,"X^{+}>=1 + X^{-}>=1");
	//hevtsum->SetDirectory(dirout);
}

ProcSkimQ::~ProcSkimQ()
{
	delete hevtsum;
	delete hists;
}

void ProcSkimQ::handle()
{
	//Info("ProcSkimQ::handle()", "");
	pass = kFALSE;
		
	hevtsum->Fill(EVT);
	if ( dH10->id[0]==ELECTRON ) {
		hevtsum->Fill(EVT_1ST_E);
		Int_t numPos=0;
		Int_t numNeg=0;
		Int_t num0=0;
		for (Int_t i=1;i<dH10->gpart;i++) {
		    if (dH10->q[i]==1)  numPos++;
		    if (dH10->q[i]==-1) numNeg++;
		    if (dH10->q[i]==0)  num0++;
		}
		Int_t binNum = 0;
		if (numPos>=2 && numNeg>=1) {	
			binNum=hevtsum->Fill(EVT_ETGT_2POS_ETGT_1NEG);
			dAna->skimq.isEVT_ETGT_2POS_ETGT_1NEG=kTRUE;
		}
		if (numPos>=2) {
			binNum=hevtsum->Fill(EVT_ETGT_2POS);
			dAna->skimq.isEVT_ETGT_2POS=kTRUE;
		}
		if (numPos>=1 && numNeg>=1) {
			binNum=hevtsum->Fill(EVT_ETGT_1POS_ETGT_1NEG);
			dAna->skimq.isEVT_ETGT_1POS_ETGT_1NEG=kTRUE;
		}
		if (binNum>0) {
			pass=kTRUE;
			EpProcessor::handle();
		}
	}
}

#endif // PROCSKIMQ_H
