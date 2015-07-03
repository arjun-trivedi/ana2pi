TCanvas *c[100];
void comp(char* f1name, char* f2name, bool phi_proj=kFALSE){

  TFile* f1=TFile::Open(f1name);
  y1=(THnSparse*)f1->Get("delast/cut/yield");
  if (phi_proj){
   int phi_bin_min=y1->GetAxis(1)->FindBin(116+1);
   int phi_bin_max=y1->GetAxis(1)->FindBin(118-1);
   y1->GetAxis(1)->SetRange(phi_bin_min,phi_bin_max);
  }
  TH1F* htheta1=(TH1F*)y1->Projection(0,"E");
  //htheta1->Draw();

  TFile* f2=TFile::Open(f2name);
  y2=(THnSparse*)f2->Get("delast/cut/yield");
  if (phi_proj){
   int phi_bin_min=y2->GetAxis(1)->FindBin(116+1);
   int phi_bin_max=y2->GetAxis(1)->FindBin(118-1);
   y2->GetAxis(1)->SetRange(phi_bin_min,phi_bin_max);
  }
  TH1F* htheta2=(TH1F*)y2->Projection(0,"E");

  TH1F* hdiv=(TH1F*)htheta1->Clone("hdiv");
  hdiv->SetTitle(TString::Format("%s/%s",f1->GetName(),f2->GetName()));
  hdiv->SetMinimum(0.5);
  hdiv->SetMaximum(2.0);
  hdiv->Divide(htheta2);
  int i=int(gRandom->Uniform(1,100));
  c[i]=new TCanvas(TString::Format("c%d",i),TString::Format("c%d",i));
  hdiv->Draw();
}
