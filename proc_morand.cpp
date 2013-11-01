#include "proc_morand.h"

ProcMorand::ProcMorand(TDirectory *td)
{
}

ProcMorand::~ProcMorand()
{
}

void ProcMorand::handle(DataOmega* d)
{
        if (d->threePi.W>1.8 && d->eKin.E1>0.8 && d->threePi.mm2ppip>0.1) {
		EpProcessor::handle(d);
	}
}
