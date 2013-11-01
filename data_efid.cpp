#include "data_efid.h"

DataEFid::DataEFid()
{
}

DataEFid::~DataEFid()
{
}

void DataEFid::Clear()
{
	fidE = kFALSE;
        sector = phi = theta = dc_xsc = dc_ysc = ech_x = ech_y = 0;
}
