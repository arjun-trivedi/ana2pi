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
	2. The total number of Reconstructed tracks is < 8, AND
	3. The charges of the Reconstructed tracks are identified as either:
		+ 2POS+1NEG or 2POS or 1POS1NEG
*******************************************************/


class ProcSkimQ : public EpProcessor
{

public:
	ProcSkimQ(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcSkimQ();
	
	void handle();
	
protected:
	/*static const Int_t NUM_EVTCUTS = 8;
	enum { EVT_NULL, EVT, EVT_1ST_NEG, EVT_LT8, EVT_1ST_E, EVT_1POS, EVT_2POS_EX, EVT_1POS1NEG_EX, EVT_2POS1NEG_EX };*/
	static const Int_t NUM_EVTCUTS=6;
	enum { EVT_NULL, EVT, EVT_1ST_E, EVT_LT8, EVT_2POS1NEG_EX, EVT_2POS_EX, EVT_1POS1NEG_EX };
};

ProcSkimQ::ProcSkimQ(TDirectory *td,DataH10* dataH10,DataAna* dataAna) 
					:EpProcessor(td, dataH10, dataAna)
{
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	//hevtsum->SetMinimum(0);
	/*hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1ST_NEG,"X^{-} 1^{st}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_LT8,"gpart < 8");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1ST_E,"e^{-} 1^{st}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1POS,"e^{-} + X^{+} + X");
	hevtsum->GetXaxis()->SetBinLabel(EVT_2POS_EX,"e^{-} + 2X^{+}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1POS1NEG_EX,"e^{-} + X^{+} + X^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_2POS1NEG_EX,"e^{-} + 2X^{+} + X^{-}");*/

	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1ST_E,"e^{-} 1^{st}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_LT8,"gpart < 8");
	hevtsum->GetXaxis()->SetBinLabel(EVT_2POS1NEG_EX,"2X^{+} + X^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_2POS_EX,"2X^{+}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1POS1NEG_EX,"X^{+} + X^{-}");
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
		//if ( dH10->gpart<8) { 
		    hevtsum->Fill(EVT_LT8);
			Int_t numPos=0;
			Int_t numNeg=0;
			Int_t num0=0;
			for (Int_t i=1;i<dH10->gpart;i++) {
			    if      (dH10->q[i]==1)  numPos++;
			    else if (dH10->q[i]==-1) numNeg++;
			    else if (dH10->q[i]==0)  num0++;
			}
			Int_t binNum = 0;
			if (numPos>=2) {//if (numPos==2 && numNeg==0) {	
				binNum=hevtsum->Fill(EVT_2POS_EX);
				dAna->skimq.isEVT_2POS_EX=kTRUE;
			}
			if (numPos>=1 && numNeg>=1) {//if (numPos==1 && numNeg==1) {
				binNum=hevtsum->Fill(EVT_1POS1NEG_EX);
				dAna->skimq.isEVT_1POS1NEG_EX=kTRUE;
			}
			if (numPos>=2 && numNeg>=1) {//if (numPos==2 && numNeg==1) {	
				binNum=hevtsum->Fill(EVT_2POS1NEG_EX);
				dAna->skimq.isEVT_2POS1NEG_EX=kTRUE;
			}
			if (binNum>0) {
				pass=kTRUE;
				EpProcessor::handle();
			}
		//} 
	}
}

#endif // PROCSKIMQ_H
