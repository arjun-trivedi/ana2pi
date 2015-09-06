#define h10looper_e1f_cxx
#include "h10looper_e1f.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void h10looper_e1f::Loop()
{
   if (fChain == 0) return;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<_nentries_to_proc;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
