#include "data_ekin.h"

DataEkin::DataEkin()
{
}

DataEkin::~DataEkin()
{
}

void DataEkin::Clear()
{
	sector = W = Q2 = nu = xb = E1 = theta1 = phi1 =theta = phi = 0;
}
