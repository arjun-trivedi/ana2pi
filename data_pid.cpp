#include "data_pid.h"

DataPid::DataPid()
{
}

DataPid::~DataPid()
{
}

void DataPid::Clear()
{
	h10IdxP = h10IdxPip = h10IdxPim =  sectorP = sectorPip = sectorPim = betaP = betaPip = betaPim = pP = pPip = pPim =  betaStrP = betaStrPip = betaStrPim = dtP = dtPip = dtPim = 0;
        P_ec_ei = P_ec_eo = P_etot = Pip_ec_ei = Pip_ec_eo = Pip_etot = Pim_ec_ei = Pim_ec_eo = Pim_etot = P_nphe = P_cc_segm = P_cc_theta = Pip_nphe = Pip_cc_segm = Pip_cc_theta = Pim_nphe = Pim_cc_segm = Pim_cc_theta = 0;
}
