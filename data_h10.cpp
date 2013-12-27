#include "data_h10.h"
#include "particle_constants.h"

using namespace ParticleConstants;

DataH10::DataH10(TString h10type)
{
	/* Determint h10type */
	is_e1f = is_e16 = is_exp = is_sim = kFALSE;
	TObjArray *h10type_tokens = h10type.Tokenize(":");
	exp   = h10type_tokens->At(0)->GetName();
	dtype = h10type_tokens->At(1)->GetName();
	if (h10type_tokens->GetEntries() == 3) skim = h10type_tokens->At(2)->GetName();
	
	if (exp.EqualTo("e1f")) {
		is_e1f = kTRUE;
		_4vE0 = E1F_4vE0;
		_4vP0 = E1F_4vP0;
	}else if (exp.EqualTo("e16")) {
		is_e16 = kTRUE;
		_4vE0 = E16_4vE0;
		_4vP0 = E16_4vP0;
	}
	else Info("DataH10::DataH10()", "Could not determine h10type.experiment!\n");
	
	if (dtype.EqualTo("exp")) is_exp = kTRUE;
	else if (dtype.EqualTo("sim")) is_sim = kTRUE;
	else Info("DataH10::DataH10()", "Could not determine h10type.dtype!\n");
	
	Info("DataH10::DataH10()", "DataH10 intitialized with following h10type: %s:%s:%s\n", exp.Data(), dtype.Data(), skim.Data());
	/* *** */
	
	run = 0;
	memset(fn, 0, 256);
	memset(bn, 0, 64);
}
void DataH10::Bind(TTree* tree)
{
	fChain = tree;
	//Bind SEB Branches
	if (_h10typ.exp=="e1f" && _h10typ.dtype=="exp") fChain->SetBranchAddress("evthel", &evthel, &b_evthel); 
	if (_h10typ.rctn=="2pi_userana" || _h10typ.rctn=="elas_userana" ||
		(_h10typ.exp=="e16" && _h10typ.dtype=="exp")) {//Bind SEB' Branches
		fChain->SetBranchAddress("id", tmpVars.id, &b_id);
		fChain->SetBranchAddress("stat", tmpVars.stat, &b_stat);
		fChain->SetBranchAddress("dc", tmpVars.dc, &b_dc);
		fChain->SetBranchAddress("cc", tmpVars.cc, &b_cc);
		fChain->SetBranchAddress("sc", tmpVars.sc, &b_sc);
		fChain->SetBranchAddress("ec", tmpVars.ec, &b_ec);
		fChain->SetBranchAddress("q", tmpVars.q, &b_q);
		fChain->SetBranchAddress("dc_sect", tmpVars.dc_sect, &b_dc_sect);
		fChain->SetBranchAddress("dc_stat", tmpVars.dc_stat, &b_dc_stat);
		fChain->SetBranchAddress("ec_stat", tmpVars.ec_stat, &b_ec_stat);
		fChain->SetBranchAddress("ec_sect", tmpVars.ec_sect, &b_ec_sect);
		fChain->SetBranchAddress("sc_sect", tmpVars.sc_sect, &b_sc_sect);
		fChain->SetBranchAddress("sc_pd", tmpVars.sc_pd, &b_sc_pd);
		fChain->SetBranchAddress("sc_stat", tmpVars.sc_stat, &b_sc_stat);
		fChain->SetBranchAddress("cc_sect", tmpVars.cc_sect, &b_cc_sect);
		fChain->SetBranchAddress("nphe", tmpVars.nphe, &b_nphe);
	}else{//Bind SEB Branches 
		fChain->SetBranchAddress("id", id, &b_id);
		fChain->SetBranchAddress("stat", stat, &b_stat);
		fChain->SetBranchAddress("dc", dc, &b_dc);
		fChain->SetBranchAddress("cc", cc, &b_cc);
		fChain->SetBranchAddress("sc", sc, &b_sc);
		fChain->SetBranchAddress("ec", ec, &b_ec);
		fChain->SetBranchAddress("q", q, &b_q);
		fChain->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
		fChain->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
		fChain->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
		fChain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
		fChain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
		fChain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
		fChain->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
		fChain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
		fChain->SetBranchAddress("nphe", nphe, &b_nphe);
	}
	fChain->SetBranchAddress("evntid", &evntid, &b_evntid);
	fChain->SetBranchAddress("npart", &npart, &b_npart);
	fChain->SetBranchAddress("q_l", &q_l, &b_q_l);
	fChain->SetBranchAddress("t_l", &t_l, &b_t_l);
	fChain->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
	fChain->SetBranchAddress("gpart", &gpart, &b_gpart);
	fChain->SetBranchAddress("p", p, &b_p);
	fChain->SetBranchAddress("m", m, &b_m);
	fChain->SetBranchAddress("b", b, &b_b);
	fChain->SetBranchAddress("cx", cx, &b_cx);
	fChain->SetBranchAddress("cy", cy, &b_cy);
	fChain->SetBranchAddress("cz", cz, &b_cz);
	fChain->SetBranchAddress("vx", vx, &b_vx);
	fChain->SetBranchAddress("vy", vy, &b_vy);
	fChain->SetBranchAddress("vz", vz, &b_vz);
	fChain->SetBranchAddress("dc_part", &dc_part, &b_dc_part);
	fChain->SetBranchAddress("dc_xsc", dc_xsc, &b_dc_xsc);
	fChain->SetBranchAddress("dc_ysc", dc_ysc, &b_dc_ysc);
	fChain->SetBranchAddress("dc_zsc", dc_zsc, &b_dc_zsc);
	fChain->SetBranchAddress("dc_cxsc", dc_cxsc, &b_dc_cxsc);
	fChain->SetBranchAddress("dc_cysc", dc_cysc, &b_dc_cysc);
	fChain->SetBranchAddress("dc_czsc", dc_czsc, &b_dc_czsc);
	fChain->SetBranchAddress("ec_part", &ec_part, &b_ec_part);
	fChain->SetBranchAddress("etot", etot, &b_etot);
	fChain->SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
	fChain->SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
	fChain->SetBranchAddress("ech_x", ech_x, &b_ech_x);
	fChain->SetBranchAddress("ech_y", ech_y, &b_ech_y);
	fChain->SetBranchAddress("ech_z", ech_z, &b_ech_z);
	fChain->SetBranchAddress("sc_part", &sc_part, &b_sc_part);
	fChain->SetBranchAddress("sc_t", sc_t, &b_sc_t);
	fChain->SetBranchAddress("sc_r", sc_r, &b_sc_r);
	fChain->SetBranchAddress("cc_part", &cc_part, &b_cc_part);
	fChain->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
	
	//MCTK Branches
	if (_h10typ.dtyp=="sim" && _h10typ.rctn="2pi") {
		fChain->SetBranchAddress("vidmvrt", &vidmvrt, &b_vidmvrt);
		fChain->SetBranchAddress("ntrmvrt", &ntrmvrt, &b_ntrmvrt);
		fChain->SetBranchAddress("xmvrt", &xmvrt, &b_xmvrt);
		fChain->SetBranchAddress("ymvrt", &ymvrt, &b_ymvrt);
		fChain->SetBranchAddress("zmvrt", &zmvrt, &b_zmvrt);
		fChain->SetBranchAddress("ch2mvrt", &ch2mvrt, &b_ch2mvrt);
		fChain->SetBranchAddress("cxxmvrt", &cxxmvrt, &b_cxxmvrt);
		fChain->SetBranchAddress("cxymvrt", &cxymvrt, &b_cxymvrt);
		fChain->SetBranchAddress("cxzmvrt", &cxzmvrt, &b_cxzmvrt);
		fChain->SetBranchAddress("cyymvrt", &cyymvrt, &b_cyymvrt);
		fChain->SetBranchAddress("cyzmvrt", &cyzmvrt, &b_cyzmvrt);
		fChain->SetBranchAddress("stamvrt", &stamvrt, &b_stamvrt);
		fChain->SetBranchAddress("mcnentr", &mcnentr, &b_mcnentr);
		fChain->SetBranchAddress("mcnpart", &mcnpart, &b_mcnpart);
		fChain->SetBranchAddress("mcst", mcst, &b_mcst);
		fChain->SetBranchAddress("mcid", mcid, &b_mcid);
		fChain->SetBranchAddress("mcpid", mcpid, &b_mcpid);
		fChain->SetBranchAddress("mctheta", mctheta, &b_mctheta);
		fChain->SetBranchAddress("mcphi", mcphi, &b_mcphi);
		fChain->SetBranchAddress("mcp", mcp, &b_mcp);
		fChain->SetBranchAddress("mcm", mcm, &b_mcm);
		fChain->SetBranchAddress("mcvx", mcvx, &b_mcvx);
		fChain->SetBranchAddress("mcvy", mcvy, &b_mcvy);
		fChain->SetBranchAddress("mcvz", mcvz, &b_mcvz);
		fChain->SetBranchAddress("mctof", mctof, &b_mctof);
	}
	//PART Branches
	if (_h10typ.dtyp=="sim" && _h10typ.rctn="elas") {
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
//	memset(lec, 0, sizeof(UChar_t) * _MAX_PARTS);
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
//	memset(dc_trk, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(dc_stat, 0, sizeof(Char_t) * _MAX_PARTS);
	memset(dc_xsc, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(dc_ysc, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(dc_zsc, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(dc_cxsc, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(dc_cysc, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(dc_czsc, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(dc_xec, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(dc_yec, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(dc_zec, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(dc_thcc, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(dc_c2, 0, sizeof(Float_t) * _MAX_PARTS);
	ec_part = 0;
	memset(ec_stat, 0, sizeof(UShort_t) * _MAX_PARTS);
	memset(ec_sect, 0, sizeof(UChar_t) * _MAX_PARTS);
//	memset(ec_whol, 0, sizeof(Int_t) * _MAX_PARTS);
//	memset(ec_inst, 0, sizeof(Int_t) * _MAX_PARTS);
//	memset(ec_oust, 0, sizeof(Int_t) * _MAX_PARTS);
	memset(etot, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(ec_ei, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(ec_eo, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(ech_x, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(ech_y, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(ech_z, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(ec_t, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(ec_r, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(ech_x, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(ech_y, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(ech_z, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(ec_m2, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(ec_m3, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(ec_m4, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(ec_c2, 0, sizeof(Float_t) * _MAX_PARTS);
	sc_part = 0;
	memset(sc_sect, 0, sizeof(UChar_t) * _MAX_PARTS);
//	memset(sc_hit, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(sc_pd, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(sc_stat, 0, sizeof(UChar_t) * _MAX_PARTS);
//	memset(edep, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(sc_t, 0, sizeof(Float_t) * _MAX_PARTS);
	memset(sc_r, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(sc_c2, 0, sizeof(Float_t) * _MAX_PARTS);
	cc_part = 0;
	memset(cc_sect, 0, sizeof(UChar_t) * _MAX_PARTS);
//	memset(cc_hit, 0, sizeof(UChar_t) * _MAX_PARTS);
	memset(cc_segm, 0, sizeof(Int_t) * _MAX_PARTS);
	memset(nphe, 0, sizeof(UShort_t) * _MAX_PARTS);
//	memset(cc_t, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(cc_r, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(cc_c2, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(lec_sect, 0, sizeof(Int_t) * _MAX_PARTS);
//	memset(lec_hit, 0, sizeof(Int_t) * _MAX_PARTS);
//	memset(lec_stat, 0, sizeof(Int_t) * _MAX_PARTS);
//	memset(lec_etot, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(lec_t, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(lec_r, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(lec_x, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(lec_y, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(lec_z, 0, sizeof(Float_t) * _MAX_PARTS);
//	memset(lec_c2, 0, sizeof(Float_t) * _MAX_PARTS);

	mcnpart = mcnentr = stamvrt = cyzmvrt = cyymvrt = cxzmvrt = cxymvrt
	= cxxmvrt = ch2mvrt = zmvrt = ymvrt = xmvrt = ntrmvrt = vidmvrt = 0;
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
    
    //tmp variables for data type compatibility
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
