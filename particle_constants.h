#ifndef PARTICLECONSTANTS_H_
#define PARTICLECONSTANTS_H_

#include <map>
#include <TMath.h>
#include <TLorentzVector.h>

using namespace std;

namespace ParticleConstants {
	
	static const Float_t SOL = 29.9792458;
	
	//particle codes, usually PDG codes, but always those used in BOS
	static const Int_t PROTON = 2212;
	static const Int_t NEUTRON = 2112;
	static const Int_t PIP = 211;
	static const Int_t PIM = -211;
	static const Int_t PI0 = 111;
	static const Int_t KP = 321;
	static const Int_t KM = -321;
	static const Int_t PHOTON = 22;
	static const Int_t ELECTRON = 11;

	//PDG particle masses in GeV/c2
	static const Double_t MASS_P = 0.93827203;
	static const Double_t MASS_N = 0.93956556;
	static const Double_t MASS_E = 0.000511;
	static const Double_t MASS_PIP = 0.13957018;
	static const Double_t MASS_PIM = 0.13957018;
	static const Double_t MASS_PI0 = 0.1349766;
	static const Double_t MASS_KP = 0.493677;
	static const Double_t MASS_KM = 0.493677;
	static const Double_t MASS_G = 0.0;
	static const Double_t MASS_OMEGA = 0.78265; //

	typedef map<int,double> partMap;
	static partMap PARTMAP;

	static Double_t GetPdgMass(Int_t pdgId) {
		if (PARTMAP.empty()) {
			PARTMAP[PROTON] = MASS_P;
			PARTMAP[NEUTRON] = MASS_N;
			PARTMAP[ELECTRON] = MASS_E;
			PARTMAP[PIP] = MASS_PIP;
			PARTMAP[PIM] = MASS_PIM;
			PARTMAP[PI0] = MASS_PI0;
			PARTMAP[KP] = MASS_KP;
			PARTMAP[KM] = MASS_KM;
			PARTMAP[PHOTON] = MASS_G;
		}
		return PARTMAP[pdgId];
	}

    static const Float_t E1F_P = 5.497;//5.497;//5.490; 
    static const TLorentzVector E1F_lvE0(0,0,E1F_P,TMath::Sqrt(E1F_P*E1F_P+MASS_E*MASS_E));
	static const TLorentzVector E1F_lvP0(0,0,0,MASS_P);

    static const Float_t E16_P = 5.754; //changed from 5.700 post parkkj conversation on mtg 02-25-13
    static const TLorentzVector E16_lvE0(0,0,E16_P,TMath::Sqrt(E16_P*E16_P+MASS_E*MASS_E));
    static const TLorentzVector E16_lvP0(0,0,0,MASS_P);
}

namespace AnalysisConstants{
	static const Int_t NTOPS    = 4;
	static const Int_t NEVTSELS = 5; //NTOPS(2pi evt)+EVTINC(ep->X)
	enum {TOP1, TOP2, TOP3, TOP4, EVTINC}; //NOTE: EVTINC not in logical progression
																		   //

	static const Int_t NSECTORS = 7;
	enum {SECTOR0, SECTOR1, SECTOR2, SECTOR3, SECTOR4, SECTOR5, SECTOR6};
        
}

#endif /* PARTICLECONSTANTS_H_ */