#ifndef DATASCPD_H
#define DATASCPD_H

#include <TROOT.h>

class DataScpd
{

public:
	DataScpd();
	virtual ~DataScpd();
	void Clear();
	Int_t sector_e;
	Int_t sector_p;
	Int_t sector_pip;
	Float_t theta_e;
	Float_t theta_p;
	Float_t theta_pip;
	Float_t p_e;
	Float_t p_p;
	Float_t p_pip;
	Int_t sc_pd_e;
	Int_t sc_pd_p;
	Int_t sc_pd_pip;

};

#endif // DATASCPD_H
