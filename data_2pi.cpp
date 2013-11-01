#include "data_2pi.h"

Data2pi::Data2pi()
{
}

Data2pi::~Data2pi()
{
}

void Data2pi::Clear()
{
	mm2ppippim = mm2ppip = mm2ppim = mm2pippim =  W = Q2 = h = 0;
	asgnmt1.M1 = asgnmt1.M2 = asgnmt1.theta = asgnmt1.phi = asgnmt1.alpha = 0;
	asgnmt2.M1 = asgnmt2.M2 = asgnmt2.theta = asgnmt2.phi = asgnmt2.alpha = 0;
	asgnmt3.M1 = asgnmt3.M2 = asgnmt3.theta = asgnmt3.phi = asgnmt3.alpha = 0;
}
