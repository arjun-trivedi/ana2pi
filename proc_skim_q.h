#ifndef PROCSKIMQ_H
#define PROCSKIMQ_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TH1.h>

using namespace ParticleConstants;

class ProcSkimQ : public EpProcessor
{

public:
	ProcSkimQ(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcSkimQ();
	
	void handle();
	
protected:
	static const Int_t NUM_EVTCUTS = 8;
	enum { EVT_NULL, EVT, EVT_1ST_NEG, EVT_LT8, EVT_1ST_E, EVT_1POS, EVT_2POS_EX, EVT_1POS1NEG_EX, EVT_2POS1NEG_EX };
};

ProcSkimQ::ProcSkimQ(TDirectory *td,DataH10* dataH10,DataAna* dataAna) 
					:EpProcessor(td, dataH10, dataAna)
{
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1ST_NEG,"X^{-} 1^{st}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_LT8,"gpart < 8");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1ST_E,"e^{-} 1^{st}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1POS,"e^{-} + X^{+} + X");
	hevtsum->GetXaxis()->SetBinLabel(EVT_2POS_EX,"e^{-} + 2X^{+}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_1POS1NEG_EX,"e^{-} + X^{+} + X^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_2POS1NEG_EX,"e^{-} + 2X^{+} + X^{-}");
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
	if ( dH10->q[0] == -1 ) {
		hevtsum->Fill(EVT_1ST_NEG);
		//if ( dH10->gpart == 3 || dH10->gpart == 4) {
		if ( dH10->gpart < 8) { //atrivedi: 071313 
		    hevtsum->Fill(EVT_LT8);
			if ( dH10->id[0] == ELECTRON ) {
				hevtsum->Fill(EVT_1ST_E);
				Int_t numPos = 0;
				Int_t numNeg = 0;
				Int_t num0 = 0;
				for (Int_t gpnum = 1; gpnum < dH10->gpart; gpnum++) {
				    if ( dH10->q[gpnum] == 1 ) numPos++;
				    else if ( dH10->q[gpnum] == -1 ) numNeg++;
				    else if ( dH10->q[gpnum] == 0 ) num0++;
				}
			    if ( numPos > 0 ) hevtsum->Fill(EVT_1POS);
				Int_t binNum = 0;
				//if ( dH10->gpart == 3 && numPos == 2 && numNeg == 0 ) {
				if ( numPos == 2 && numNeg == 0 ) {//atrivedi: 071513 	
					binNum = hevtsum->Fill(EVT_2POS_EX);
					dAna->skimq.isEVT_2POS_EX = kTRUE;
				}
				//if ( dH10->gpart == 3 && numPos == 1 && numNeg == 1 ) {
				if ( numPos == 1 && numNeg == 1 ) {//atrivedi: 071513
					binNum = hevtsum->Fill(EVT_1POS1NEG_EX);
					dAna->skimq.isEVT_1POS1NEG_EX = kTRUE;
				}
				//if ( dH10->gpart == 4 && numPos == 2 && numNeg == 1 ) {
				if ( numPos == 2 && numNeg == 1 ) {//atrivedi: 071513	
					binNum = hevtsum->Fill(EVT_2POS1NEG_EX);
					dAna->skimq.isEVT_2POS1NEG_EX = kTRUE;
				}
				if ( binNum > 0 ) {
					pass = kTRUE;
					EpProcessor::handle();
				}
			}
		} 
	}
}

#endif // PROCSKIMQ_H
