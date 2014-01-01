#define H10Looper_cxx
#include "h10looper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void H10Looper::Loop(Long64_t nentries)
{
   if (fChain == 0) return;

   //Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      //Clear event data
      dH10->Clear();
      dAna->Clear();

      //GetEntry
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //If needed, Reconcile dH10
      if (dH10->rctn=="2pi_userana" || 
          dH10->rctn=="elast_userana" ||
          (dH10->expt=="e16" && dH10->dtyp=="exp")) {
         dH10->Reconcile();
      }
      // if (Cut(ientry) < 0) continue;
      
      Info("Loop", "Processing entry# %d\n",jentry);
      //Reset proc_chain
      //Call proc_chain
      proc_chain->handle();
   }
}
