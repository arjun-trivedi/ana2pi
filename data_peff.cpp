#include "data_peff.h"

DataPEff::DataPEff()
{
}

DataPEff::~DataPEff()
{
}

void DataPEff::Clear()
{
	sector_p=sector_pip=sector_pim-1;
	theta_p=theta_pip=theta_pim=0;
	p_p=p_pip=p_pim=0;
}
