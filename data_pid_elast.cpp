#include "data_pid_elast.h"

DataPidElast::DataPidElast()
{
}

DataPidElast::~DataPidElast()
{
}

void DataPidElast::Clear()
{
	h10IdxP = sectorP = betaP = pP = betaStrP = dtP = 0;
	P_ec_ei = P_ec_eo = P_etot = P_nphe = P_cc_segm = P_cc_theta  = 0;
}
