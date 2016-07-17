#include "data_scpd.h"

DataScpd::DataScpd()
{
}

DataScpd::~DataScpd()
{
}

void DataScpd::Clear()
{
	sector_e=sector_p=sector_pip=-1;
	theta_e=theta_p=theta_pip=0;
	p_e=p_p=p_pip=0;
	sc_pd_e=sc_pd_p=sc_pd_pip=0;
}
