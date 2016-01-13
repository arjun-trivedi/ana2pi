#include "data_mom.h"
/*
[01-10-16]
Added variables for pcorr:{pip,p}
*/
DataMom::DataMom()
{
}

DataMom::~DataMom()
{
}

void DataMom::Clear()
{
	p = sector = dcx = dcy = dcz = dp = 0;
	p_pip = sector_pip = dcx_pip = dcy_pip = dcz_pip = dp_pip = 0;
	p_p = sector_p = dcx_p = dcy_p = dcz_p = dp_p = 0;
}
