/**********************************/
/* Class that holds data structures
   to which the Branches of h10-Trees are Bound.
   For more on Branches/Banks and Binding read:
   at-h10/formats/README
*/
#ifndef DATAH10
#define DATAH10

#include <TROOT.h>
#include <TChain.h>
#include <TLorentzVector.h>

/*struct H10Typ{
	TString expt;
	TString dtyp;
	TString rctn;
};*/

class DataH10
{

public:

	DataH10(TString h10type);
	~DataH10();
	void Bind(TTree* tree);
	void Reconcile();
	void Clear();

	TString expt,dtyp,rctn;
	TTree* h10chain;
			
	Int_t run;
	Char_t fn[256];
	Char_t bn[64];
		
	TLorentzVector lvE0;
	TLorentzVector lvP0;

	static const int _MAX_PARTS = 40;
	
	/* *** Data Structures for Binding Branches *** */
	// SEB Banks
	UInt_t evntid;
	Char_t evthel;
	UChar_t npart;
	Float_t q_l;
	Float_t t_l;
	Float_t tr_time;
	Int_t gpart;
	Short_t id[_MAX_PARTS]; //[gpart]
	Char_t stat[_MAX_PARTS]; //[gpart]
	UChar_t dc[_MAX_PARTS]; //[gpart]
	UChar_t cc[_MAX_PARTS]; //[gpart]
	UChar_t sc[_MAX_PARTS]; //[gpart]
	UChar_t ec[_MAX_PARTS]; //[gpart]
	Float_t p[_MAX_PARTS]; //[gpart]
	Float_t m[_MAX_PARTS]; //[gpart]
	Char_t q[_MAX_PARTS]; //[gpart]
	Float_t b[_MAX_PARTS]; //[gpart]
	Float_t cx[_MAX_PARTS]; //[gpart]
	Float_t cy[_MAX_PARTS]; //[gpart]
	Float_t cz[_MAX_PARTS]; //[gpart]
	Float_t vx[_MAX_PARTS]; //[gpart]
	Float_t vy[_MAX_PARTS]; //[gpart]
	Float_t vz[_MAX_PARTS]; //[gpart]
	Int_t dc_part;
	UChar_t dc_sect[_MAX_PARTS]; //[dc_part]
	Char_t dc_stat[_MAX_PARTS]; //[dc_part]
	Float_t dc_xsc[_MAX_PARTS]; //[dc_part]
	Float_t dc_ysc[_MAX_PARTS]; //[dc_part]
	Float_t dc_zsc[_MAX_PARTS]; //[dc_part]
	Float_t dc_cxsc[_MAX_PARTS]; //[dc_part]
	Float_t dc_cysc[_MAX_PARTS]; //[dc_part]
	Float_t dc_czsc[_MAX_PARTS]; //[dc_part]
	Int_t ec_part;
	UShort_t ec_stat[_MAX_PARTS]; //[ec_part]
	UChar_t ec_sect[_MAX_PARTS]; //[ec_part]
	Float_t etot[_MAX_PARTS]; //[ec_part]
	Float_t ec_ei[_MAX_PARTS]; //[ec_part]
	Float_t ec_eo[_MAX_PARTS]; //[ec_part]
	Float_t ech_x[_MAX_PARTS]; //[ec_part]
	Float_t ech_y[_MAX_PARTS]; //[ec_part]
	Float_t ech_z[_MAX_PARTS]; //[ec_part]
	Int_t sc_part;
	UChar_t sc_sect[_MAX_PARTS]; //[sc_part]
	UChar_t sc_pd[_MAX_PARTS]; //[sc_part]
	UChar_t sc_stat[_MAX_PARTS]; //[sc_part]
	Float_t sc_t[_MAX_PARTS]; //[sc_part]
	Float_t sc_r[_MAX_PARTS]; //[sc_part]
	Int_t cc_part;
	UChar_t cc_sect[_MAX_PARTS]; //[cc_part]
	Int_t cc_segm[_MAX_PARTS]; //[cc_part]
	UShort_t nphe[_MAX_PARTS]; //[cc_part]
	
	//MCTK Banks
	Int_t mcnentr;
	UChar_t mcnpart;
	Int_t mcst[_MAX_PARTS]; //[mcnentr]
	Int_t mcid[_MAX_PARTS]; //[mcnentr]
	Int_t mcpid[_MAX_PARTS]; //[mcnentr]
	Float_t mctheta[_MAX_PARTS]; //[mcnentr]
	Float_t mcphi[_MAX_PARTS]; //[mcnentr]
	Float_t mcp[_MAX_PARTS]; //[mcnentr]
	Float_t mcm[_MAX_PARTS]; //[mcnentr]
	Float_t mcvx[_MAX_PARTS]; //[mcnentr]
	Float_t mcvy[_MAX_PARTS]; //[mcnentr]
	Float_t mcvz[_MAX_PARTS]; //[mcnentr]
	Float_t mctof[_MAX_PARTS]; //[mcnentr]

	//PART Banks
	Int_t nprt;
	Int_t pidpart[_MAX_PARTS];   //[nprt]
	Float_t xpart[_MAX_PARTS];   //[nprt]
	Float_t ypart[_MAX_PARTS];   //[nprt]
	Float_t zpart[_MAX_PARTS];   //[nprt]
	Float_t epart[_MAX_PARTS];   //[nprt]
	Float_t pxpart[_MAX_PARTS];   //[nprt]
	Float_t pypart[_MAX_PARTS];   //[nprt]
	Float_t pzpart[_MAX_PARTS];   //[nprt]
	Float_t qpart[_MAX_PARTS];   //[nprt]
	Int_t flagspart[_MAX_PARTS];   //[nprt]
	
	//for some Branches of SEB' Banks that are different
	struct tmp{
		Int_t id[_MAX_PARTS]; //[gpart]
		Int_t stat[_MAX_PARTS]; //[gpart]
		Int_t dc[_MAX_PARTS]; //[gpart]
		Int_t cc[_MAX_PARTS]; //[gpart]
		Int_t sc[_MAX_PARTS]; //[gpart]
		Int_t ec[_MAX_PARTS]; //[gpart]
		Int_t q[_MAX_PARTS]; //[gpart]
		Int_t dc_sect[_MAX_PARTS]; //[dc_part]
		Int_t dc_stat[_MAX_PARTS]; //[dc_part]
		Int_t ec_stat[_MAX_PARTS]; //[ec_part]
		Int_t ec_sect[_MAX_PARTS]; //[ec_part]
		Int_t sc_sect[_MAX_PARTS]; //[sc_part]
		Int_t sc_pd[_MAX_PARTS]; //[sc_part]
		Int_t sc_stat[_MAX_PARTS]; //[sc_part]
		Int_t cc_sect[_MAX_PARTS]; //[cc_part]
		Int_t nphe[_MAX_PARTS]; //[cc_part]
	} tmpVars;
	
		
	/* *** Branches in h10 Trees *** */
	//SEB (=SEB') Branches
	TBranch *b_evntid; //!
	TBranch *b_evthel;
	TBranch *b_npart; //!
	TBranch *b_q_l; //!
	TBranch *b_t_l; //!
	TBranch *b_tr_time; //!
	TBranch *b_gpart; //!
	TBranch *b_id; //!
	TBranch *b_stat; //!
	TBranch *b_dc; //!
	TBranch *b_cc; //!
	TBranch *b_sc; //!
	TBranch *b_ec; //!
	TBranch *b_p; //!
	TBranch *b_m; //!
	TBranch *b_q; //!
	TBranch *b_b; //!
	TBranch *b_cx; //!
	TBranch *b_cy; //!
	TBranch *b_cz; //!
	TBranch *b_vx; //!
	TBranch *b_vy; //!
	TBranch *b_vz; //!
	TBranch *b_dc_part; //!
	TBranch *b_dc_sect; //!
	TBranch *b_dc_stat; //!
	TBranch *b_dc_xsc; //!
	TBranch *b_dc_ysc; //!
	TBranch *b_dc_zsc; //!
	TBranch *b_dc_cxsc; //!
	TBranch *b_dc_cysc; //!
	TBranch *b_dc_czsc; //!
	TBranch *b_ec_part; //!
	TBranch *b_ec_stat; //!
	TBranch *b_ec_sect; //!
	TBranch *b_etot; //!
	TBranch *b_ec_ei; //!
	TBranch *b_ec_eo; //!
	TBranch *b_ech_x; //!
	TBranch *b_ech_y; //!
	TBranch *b_ech_z; //!
	TBranch *b_sc_part; //!
	TBranch *b_sc_sect; //!
	TBranch *b_sc_pd; //!
	TBranch *b_sc_stat; //!
	TBranch *b_sc_t; //!
	TBranch *b_sc_r; //!
	TBranch *b_cc_part; //!
	TBranch *b_cc_sect; //!
	TBranch *b_cc_segm; //!
	TBranch *b_nphe; //!

	//MCTK Branches 
	TBranch *b_mcnentr; //!
	TBranch *b_mcnpart; //!
	TBranch *b_mcst; //!
	TBranch *b_mcid; //!
	TBranch *b_mcpid; //!
	TBranch *b_mctheta; //!
	TBranch *b_mcphi; //!
	TBranch *b_mcp; //!
	TBranch *b_mcm; //!
	TBranch *b_mcvx; //!
	TBranch *b_mcvy; //!
	TBranch *b_mcvz; //!
	TBranch *b_mctof; //!

	//PART Branches
	TBranch *b_nprt;   //!
	TBranch *b_pidpart;   //!
	TBranch *b_xpart;   //!
	TBranch *b_ypart;   //!
	TBranch *b_zpart;   //!
	TBranch *b_epart;   //!
	TBranch *b_pxpart;   //!
	TBranch *b_pypart;   //!
	TBranch *b_pzpart;   //!
	TBranch *b_qpart;   //!
	TBranch *b_flagspart;   //!

};

#endif // DATAH10
