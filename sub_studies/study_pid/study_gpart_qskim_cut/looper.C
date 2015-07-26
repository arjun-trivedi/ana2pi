#define looper_cxx
#include "looper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void looper::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L looper.C
//      Root > looper t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   int NCUTS=5;
   enum {EVT_NULL,P1_PIP1,P1_PIP2,P2_PIP1,P2_PIP2,OTHER};
   TH1F* h[4];
   h[0]=new TH1F("h_ex_2pos","npos=2",NCUTS,0.5,NCUTS+0.5);
   h[1]=new TH1F("h_gt_2pos","npos>2",NCUTS,0.5,NCUTS+0.5);
   h[2]=new TH1F("h_ex_2pos_gpart_lt_8","npos=2 & gpart<8",NCUTS,0.5,NCUTS+0.5);
   h[3]=new TH1F("h_gt_2pos_gpart_lt_8","npos>2 & gpart<8",NCUTS,0.5,NCUTS+0.5);
  
   for (int i=0;i<4;i++){ 
      h[i]->GetXaxis()->SetBinLabel(P1_PIP1,"1p_1pip");
      h[i]->GetXaxis()->SetBinLabel(P1_PIP2,"1p_2pip");
      h[i]->GetXaxis()->SetBinLabel(P2_PIP1,"2p_1pip");
      h[i]->GetXaxis()->SetBinLabel(P2_PIP2,"2p_2pip");
      h[i]->GetXaxis()->SetBinLabel(OTHER,"other");
   }

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //cout <<"Processing entry "<<jentry<<endl;
      if( (np==1) && (npip==1) ){
         if (npos==2) h[0]->Fill(P1_PIP1);
         if (npos>2)  h[1]->Fill(P1_PIP1);
         if (npos==2 && gpart<8) h[2]->Fill(P1_PIP1);
         if (npos>2 && gpart<8)  h[3]->Fill(P1_PIP1);
      }else if ( (np==1) && (npip==2) ){
         if (npos==2) h[0]->Fill(P1_PIP2);
         if (npos>2) h[1]->Fill(P1_PIP2);
         if (npos==2 && gpart<8) h[2]->Fill(P1_PIP2);
         if (npos>2 && gpart<8)  h[3]->Fill(P1_PIP2);
      }else if ( (np==2) && (npip==1) ){
         if (npos==2) h[0]->Fill(P2_PIP1);
         if (npos>2) h[1]->Fill(P2_PIP1);
         if (npos==2 && gpart<8) h[2]->Fill(P2_PIP1);
         if (npos>2 && gpart<8)  h[3]->Fill(P2_PIP1);
      }else if ( (np==2) && (npip==2) ){
         if (npos==2) h[0]->Fill(P2_PIP2);
         if (npos>2) h[1]->Fill(P2_PIP2);
         if (npos==2 && gpart<8) h[2]->Fill(P2_PIP2);
         if (npos>2 && gpart<8)  h[3]->Fill(P2_PIP2);
      }else{
         if (npos==2) h[0]->Fill(OTHER);
         if (npos>2) h[1]->Fill(OTHER);
         if (npos==2 && gpart<8) h[2]->Fill(OTHER);
         if (npos>2 && gpart<8)  h[3]->Fill(OTHER);
      }
   }
   TCanvas *c[4];
   for (int i=0;i<4;i++){
     c[i]=new TCanvas(TString::Format("c%d",i+1),TString::Format("c%d",i+1));
     h[i]->Draw();
   }
}
