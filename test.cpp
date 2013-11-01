#include "eid.h"
#include "TCanvas.h"
void test() {
  Eid* m = new Eid();
  TObjArray* objs = m->getFuncsCuts();

  TCanvas *c_func = new TCanvas("func", "func");
  c_func->Divide(3,2);

  for (Int_t iSector=0; iSector<6; iSector++){
    c_func->cd(iSector+1);

    Int_t j = iSector*3;
    TF1* o0 = (TF1*)objs->At(j);
    o0->SetMaximum(0.5);
    o0->SetMinimum(0.);
    o0->Draw();
    objs->At(j+1)->Draw("same");
    objs->At(j+2)->Draw("same");
  }
}
