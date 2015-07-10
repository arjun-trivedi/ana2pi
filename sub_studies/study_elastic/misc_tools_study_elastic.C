TCanvas *c[100];
void comp(char* f1name, char* f2name, bool phi_proj=kFALSE,bool dbg=kFALSE,bool use_fullrun=kFALSE){

  TFile* f1=TFile::Open(f1name);
  THnSparse* y1;
  if (use_fullrun) y1=(THnSparse*)f1->Get("delast/yield");
  else             y1=(THnSparse*)f1->Get("delast/cut/yield");
  if (phi_proj){
   int phi_bin_min=y1->GetAxis(1)->FindBin(354+1);
   int phi_bin_max=y1->GetAxis(1)->FindBin(356-1);
   y1->GetAxis(1)->SetRange(phi_bin_min,phi_bin_max);
  }
  TH1F* htheta1=(TH1F*)y1->Projection(0,"E");
  if (dbg){
    htheta1->SetTitle(TString::Format("htheta_%s",f1->GetName()));
    int j=int(gRandom->Uniform(1,100));
    c[j]=new TCanvas(TString::Format("c%d",j),TString::Format("c%d",j));
    c[j]->Divide(1,2);
    c[j]->cd(1);
    htheta1->Draw();
  }

  TFile* f2=TFile::Open(f2name);
  THnSparse* y2;
  if (use_fullrun) y2=(THnSparse*)f2->Get("delast/yield");
  else             y2=(THnSparse*)f2->Get("delast/cut/yield");
  if (phi_proj){
   int phi_bin_min=y2->GetAxis(1)->FindBin(354+1);
   int phi_bin_max=y2->GetAxis(1)->FindBin(356-1);
   y2->GetAxis(1)->SetRange(phi_bin_min,phi_bin_max);
  }
  TH1F* htheta2=(TH1F*)y2->Projection(0,"E");
  if (dbg){
    htheta2->SetTitle(TString::Format("htheta_%s",f2->GetName()));
    c[j]->cd(2);
    htheta2->Draw();
  }

  TH1F* hdiv=(TH1F*)htheta1->Clone("hdiv");
  hdiv->SetTitle(TString::Format("%s/%s",f1->GetName(),f2->GetName()));
  hdiv->SetMinimum(0.5);
  hdiv->SetMaximum(2.0);
  hdiv->Divide(htheta2);
  int i=int(gRandom->Uniform(1,100));
  c[i]=new TCanvas(TString::Format("c%d",i),TString::Format("c%d",i));
  hdiv->Draw();
}
