#include "data_elastic.h"

DataElastic::DataElastic()
{
}

DataElastic::~DataElastic()
{
}

void DataElastic::Clear()
{
	gpart=ne=np=Q2=W=MMp=0;//p[0]=p[1]=px[0]=px[1]=py[0]=py[1]=pz[0]=pz[1]=0;
	lvE.SetXYZT(0,0,0,0);
	lvP.SetXYZT(0,0,0,0);
}
