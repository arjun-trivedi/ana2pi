#include "data_h10.h"
#include "particle_constants.h"

using namespace ParticleConstants;

DataH10::DataH10(TString h10type)
{
	//Determint h10typ
	TObjArray *h10type_tokens = h10type.Tokenize(":");
	h10typ.exp  = h10type_tokens->At(0)->GetName();
	h10typ.dtyp = h10type_tokens->At(1)->GetName();
	h10typ.rctn = h10type_tokens->At(2)->GetName();

	if (h10typ.exp=="e1f") {
		lvE0 = E1F_lvE0;
		lvP0 = E1F_lvP0;
	}else if (h10typ.exp=="e16") {
		lvE0 = E16_lvE0;
		lvP0 = E16_lvP0;
	}
	else Info("DataH10::DataH10()", "Could not determine h10typ.exp!\n");
		
	run = 0;
	memset(fn, 0, 256);
	memset(bn, 0, 64);
}
void DataH10::Bind(TTree* tree)
{
	h10chain = tree;
	//Bind SEB Branches
	if (h10typ.exp=="e1f" && h10typ.dtyp=="exp") h10chain->SetBranchAddress("evthel", &evthel, &b_evthel); 
	if (h10typ.rctn=="2pi_userana" || h10typ.rctn=="elas_userana" ||
		(h10typ.exp=="e16" && h10typ.dtyp=="exp")) {//Bind SEB' Branches
		h10chain->SetBranchAddress("id", tmpVars.id, &b_id);
		h10chain->SetBranchAddress("stat", tmpVars.stat, &b_stat);
		h10chain->SetBranchAddress("dc", tmpVars.dc, &b_dc);
		h10chain->SetBranchAddress("cc", tmpVars.cc, &b_cc);
		h10chain->SetBranchAddress("sc", tmpVars.sc, &b_sc);
		h10chain->SetBranchAddress("ec", tmpVars.ec, &b_ec);
		h10chain->SetBranchAddress("q", tmpVars.q, &b_q);
		h10chain->SetBranchAddress("dc_sect", tmpVars.dc_sect, &b_dc_sect);
		h10chain->SetBranchAddress("dc_stat", tmpVars.dc_stat, &b_dc_stat);
		h10chain->SetBranchAddress("ec_stat", tmpVars.ec_stat, &b_ec_stat);
		h10chain->SetBranchAddress("ec_sect", tmpVars.ec_sect, &b_ec_sect);
		h10chain->SetBranchAddress("sc_sect", tmpVars.sc_sect, &b_sc_sect);
		h10chain->SetBranchAddress("sc_pd", tmpVars.sc_pd, &b_sc_pd);
		h10chain->SetBranchAddress("sc_stat", tmpVars.sc_stat, &b_sc_stat);
		h10chain->SetBranchAddress("cc_sect", tmpVars.cc_sect, &b_cc_sect);
		h10chain->SetBranchAddress("nphe", tmpVars.nphe, &b_nphe);
	}else{//Bind SEB Branches 
		h10chain->SetBranchAddress("id", id, &b_id);
		h10chain->SetBranchAddress("stat", stat, &b_stat);
		h10chain->SetBranchAddress("dc", dc, &b_dc);
		h10chain->SetBranchAddress("cc", cc, &b_cc);
		h10chain->SetBranchAddress("sc", sc, &b_sc);
		h10chain->SetBranchAddress("ec", ec, &b_ec);
		h10chain->SetBranchAddress("q", q, &b_q);
		h10chain->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
		h10chain->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
		h10chain->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
		h10chain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
		h10chain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
		h10chain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
		h10chain->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
		h10chain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
		h10chain->SetBranchAddress("nphe", nphe, &b_nphe);
	}
	h10chain->SetBranchAddress("evntid", &evntid, &b_evntid);
	h10chain->SetBranchAddress("npart", &npart, &b_npart);
	h10chain->SetBranchAddress("q_l", &q_l, &b_q_l);
	h10chain->SetBranchAddress("t_l", &t_l, &b_t_l);
	h10chain->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
	h10chain->SetBranchAddress("gpart", &gpart, &b_gpart);
	h10chain->SetBranchAddress("p", p, &b_p);
	h10chain->SetBranchAddress("m", m, &b_m);
	h10chain->SetBranchAddress("b", b, &b_b);
	h10chain->SetBranchAddress("cx", cx, &b_cx);
	h10chain->SetBranchAddress("cy", cy, &b_cy);
	h10chain->SetBranchAddress("cz", cz, &b_cz);
	h10chain->SetBranchAddress("vx", vx, &b_vx);
	h10chain->SetBranchAddress("vy", vy, &b_vy);
	h10chain->SetBranchAddress("vz", vz, &b_vz);
	h10chain->SetBranchAddress("dc_part", &dc_part, &b_dc_part);
	h10chain->SetBranchAddress("dc_xsc", dc_xsc, &b_dc_xsc);
	h10chain->SetBranchAddress("dc_ysc", dc_ysc, &b_dc_ysc);
	h10chain->SetBranchAddress("dc_zsc", dc_zsc, &b_dc_zsc);
	h10chain->SetBranchAddress("dc_cxsc", dc_cxsc, &b_dc_cxsc);
	h10chain->SetBranchAddress("dc_cysc", dc_cysc, &b_dc_cysc);
	h10chain->SetBranchAddress("dc_czsc", dc_czsc, &b_dc_czsc);
	h10chain->SetBranchAddress("ec_part", &ec_part, &b_ec_part);
	h10chain->SetBranchAddress("etot", etot, &b_etot);
	h10chain->SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
	h10chain->SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
	h10chain->SetBranchAddress("ech_x", ech_x, &b_ech_x);
	h10chain->SetBranchAddress("ech_y", ech_y, &b_ech_y);
	h10chain->SetBranchAddress("ech_z", ech_z, &b_ech_z);
	h10chain->SetBranchAddress("sc_part", &sc_part, &b_sc_part);
	h10chain->SetBranchAddress("sc_t", sc_t, &b_sc_t);
	h10chain->SetBranchAddress("sc_r", sc_r, &b_sc_r);
	h10chain->SetBranchAddress("cc_part", &cc_part, &b_cc_part);
	h10chain->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
	
	//MCTK Branches
	if (h10typ.dtyp=="sim" && h10typ.rctn=="2pi") {
		h10chain->SetBranchAddress("mcnentr", &mcnentr, &b_mcnentr);
		h10chain->SetBranchAddress("mcnpart", &mcnpart, &b_mcnpart);
		h10chain->SetBranchAddress("mcst", mcst, &b_mcst);
		h10chain->SetBranchAddress("mcid", mcid, &b_mcid);
		h10chain->SetBranchAddress("mcpid", mcpid, &b_mcpid);
		h10chain->SetBranchAddress("mctheta", mctheta, &b_mctheta);
		h10chain->SetBranchAddress("mcphi", mcphi, &b_mcphi);
		h10chain->SetBranchAddress("mcp", mcp, &b_mcp);
		h10chain->SetBranchAddress("mcm", mcm, &b_mcm);
		h10chain->SetBranchAddress("mcvx", mcvx, &b_mcvx);
		h10chain->SetBranchAddress("mcvy", mcvy, &b_mcvy);
		h10chain->SetBranchAddress("mcvz", mcvz, &b_mcvz);
		h10chain->SetBranchAddress("mctof", mctof, &b_mctof);
	}
	//PART Branches
	if (h10typ.dtyp=="sim" && h10typ.rctn=="elas") {
	}
}

void DataH10::Reconcile() {
	for (int i = 0; i < gpart; i++) {
		id[i] = tmpVars.id[i];
		stat[i] = tmpVars.stat[i];
		dc[i] = tmpVars.dc[i];
		cc[i] = tmpVars.cc[i];
		sc[i] = tmpVars.sc[i];
		ec[i] = tmpVars.ec[i];
		q[i] = tmpVars.q[i];
	}
	for (int i = 0; i < dc_part; i++) {
		dc_sect[i] = tmpVars.dc_sect[i];
		dc_stat[i] = tmpVars.dc_stat[i];
	}
	for (int i = 0; i < ec_part; i++) {
		ec_stat[i] = tmpVars.ec_stat[i];
		ec_sect[i] = tmpVars.ec_sect[i];
	}
	for (int i = 0; i < sc_part; i++) {
		sc_stat[i] = tmpVars.sc_stat[i];
		sc_sect[i] = tmpVars.sc_sect[i];
		sc_pd[i] = tmpVars.sc_pd[i];
	}
	for (int i = 0; i < cc_part; i++) {
		cc_sect[i] = tmpVars.cc_sect[i];
		nphe[i] = tmpVars.nphe[i];
	}
}

DataH10::~DataH10()
{
}

void DataH10::Clear()
{	
	//SEB
	evntid = 0;
	evthel = 0;
	npart = 0;;
	q_l = 0;
	t_l = 0;
	tr_time = 0;
	gpart = 0;
	memset(id, 0, sizeof(Short_t) * _MAX_PARTS);
	memset(stat, 0, sizeof(Char_t) * _MAX_PARTS);
	memset(dc, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(cc, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(sc, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(ec, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(p, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(m, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(q, 0, sizeof(Char_t) * _MAX_PARTS);
	memset(b, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(cx, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(cy, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(cz, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(vx, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(vy, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(vz, 0, sizeof(Float_t) * _MAX_PARTS);
	dc_part = 0;
	memset(dc_sect, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(dc_stat, 0, sizeof(Char_t) * _MAX_PARTS);
	memset(dc_xsc, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(dc_ysc, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(dc_zsc, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(dc_cxsc, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(dc_cysc, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(dc_czsc, 0, sizeof(Float_t) * _MAX_PARTS);
	ec_part = 0;
	memset(ec_stat, 0, sizeof(UShort_t) * _MAX_PARTS);
	memset(ec_sect, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(etot, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(ec_ei, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(ec_eo, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(ech_x, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(ech_y, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(ech_z, 0, sizeof(Float_t) * _MAX_PARTS);
	sc_part = 0;
	memset(sc_sect, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(sc_pd, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(sc_stat, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(sc_t, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(sc_r, 0, sizeof(Float_t) * _MAX_PARTS);
	cc_part = 0;
	memset(cc_sect, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(cc_segm, 0, sizeof(Int_t) * _MAX_PARTS);
	memset(nphe, 0, sizeof(UShort_t) * _MAX_PARTS);

	//MCTK
	mcnentr = mcnpart = 0;
	memset(mcid, 0, sizeof(Int_t) * _MAX_PARTS);
	memset(mcst, 0, sizeof(Int_t) * _MAX_PARTS);
	memset(mcpid, 0, sizeof(Int_t) * _MAX_PARTS);
	memset(mctheta, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(mcphi, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(mcp, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(mcm, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(mcvx, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(mcvy, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(mcvz, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(mctof, 0, sizeof(Float_t) * _MAX_PARTS);
    
    //SEB'
	memset(tmpVars.id, 0, sizeof(Int_t) * _MAX_PARTS);
	memset(tmpVars.stat, 0, sizeof(Int_t) * _MAX_PARTS);
	memset(tmpVars.dc, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.cc, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.sc, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.ec, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.q, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.dc_sect, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.dc_stat, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.ec_stat, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.ec_sect, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.sc_sect, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.sc_pd, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.sc_stat, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.cc_sect, 0, sizeof(Int_t) * _MAX_PARTS); 
	memset(tmpVars.nphe, 0, sizeof(Int_t) * _MAX_PARTS); 
}
