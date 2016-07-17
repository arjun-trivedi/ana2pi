#include "data_eff.h"

DataEff::DataEff()
{
}

DataEff::~DataEff()
{
}

void DataEff::Clear()
{
	sector_e=sector_p=sector_pip=-1;
	theta_e=theta_p=theta_pip=0;
	p_e=p_p=p_pip=0;
}
