#include "data_2pievt.h"

Data2piEvt::Data2piEvt()
{
}

Data2piEvt::~Data2piEvt()
{
}

void Data2piEvt::Clear()
{
	mm2ppippim = mm2ppip = mm2ppim = mm2pippim =  W = Q2 = h = 0;
	varset1.M1 = varset1.M2 = varset1.theta = varset1.phi = varset1.alpha = 0;
	varset2.M1 = varset2.M2 = varset2.theta = varset2.phi = varset2.alpha = 0;
	varset3.M1 = varset3.M2 = varset3.theta = varset3.phi = varset3.alpha = 0;
}
