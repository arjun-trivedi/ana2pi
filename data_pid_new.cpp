#include "data_pid_new.h"

DataPidNew::DataPidNew()
{
}

DataPidNew::~DataPidNew()
{
}

void DataPidNew::Clear()
{
	h10IdxP=h10IdxPip=h10IdxPim=0;
	l_e=t_e=t_off=0;
	ntrk=0;
	for (int i=0;i<kMaxTrack;i++){
		q[i]=dc[i]=sc[i]=p[i]=l[i]=t[i]=b[i]=sector[i]=id[i]=h10_idx[i]=0;
		b_p[i]=b_pip[i]=b_pim[i]=0;
		dt_p[i]=dt_pip[i]=dt_pim[i]=0;
	}
}