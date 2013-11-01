#include "data_h10skim.h"
#include "particle_constants.h"

using namespace ParticleConstants;

DataH10skim::DataH10skim() : DataH10()
{
	
}
void DataH10skim::Bind(TTree* fChain)
{
	//Info("DataH10skim::Bind()", "binding branches\n");
	if (is_e1f) fChain->SetBranchAddress("evthel", &evthel, &b_evthel);
	fChain->SetBranchAddress("id", id, &b_id);
	fChain->SetBranchAddress("gpart", &gpart, &b_gpart);
	fChain->SetBranchAddress("p", p, &b_p);
	fChain->SetBranchAddress("cx", cx, &b_cx);
	fChain->SetBranchAddress("cy", cy, &b_cy);
	fChain->SetBranchAddress("cz", cz, &b_cz);
}

DataH10skim::~DataH10skim()
{
}

void DataH10skim::Clear()
{	
	evthel = 0;
    gpart = 0;
	memset(id, 0, sizeof(Short_t) * _MAX_PARTS);
	memset(p, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(cx, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(cy, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(cz, 0, sizeof(Float_t) * _MAX_PARTS);
}
