#ifndef PARTICLECONSTANTS_H_
#define PARTICLECONSTANTS_H_

#include <map>
#include <TMath.h>
#include <TLorentzVector.h>

using namespace std;

namespace ParticleConstants {

	static const Float_t E1F_E0 = 5.497; //e1f
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
	
	static const TLorentzVector _4vE0(0,0,E1F_E0,TMath::Sqrt(E1F_E0*E1F_E0+MASS_E*MASS_E));
	static const TLorentzVector _4vP0(0,0,0,MASS_P);
}

#endif /* PARTICLECONSTANTS_H_ */
