#include "data_top.h"

DataTop::DataTop()
{
}

DataTop::~DataTop()
{
}

void DataTop::Clear()
{
	mm2ppippim = mmppippim  = mm2ppip = mmppip= mm2ppim = mmppim = mm2pippim =  mmpippim = Q2 = W = h = 0;
	varset1.M1 = varset1.M2 = varset1.theta = varset1.phi = varset1.alpha = 0;
	varset2.M1 = varset2.M2 = varset2.theta = varset2.phi = varset2.alpha = 0;
	varset3.M1 = varset3.M2 = varset3.theta = varset3.phi = varset3.alpha = 0;
}
