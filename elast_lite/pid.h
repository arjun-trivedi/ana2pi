#ifndef PID_H
#define PID_H

#include <fstream>
#include <string>
#include <TF1.h>

using namespace std;

class Pid {

public:
	Pid(TString dtyp);
	~Pid();
	
	bool is_proton(float dt,float p);
	bool is_pip(float dt,float p);
	bool is_pim(float dt,float p);

	void plot_dtVp_cuts();
private:
	TString _dtyp;
	static const int _NPRTCLS=3;
	enum {_P,_PIP,_PIM};
	TString* _prtcl_name;

	TF1** _fcut_l;
	TF1** _fcut_h;
};

#endif // PID_H
