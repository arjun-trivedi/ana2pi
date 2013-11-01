Double_t calcError(THnSparse*  h1, THnSparse* h2){

  //Bool_t wantErrors=h1->GetCalculateErrors();
  //if (!h1->GetCalculateErrors() && h2->GetCalculateErrors())
  //    wantErrors=kTRUE;

  Int_t coord[8];
  Double_t err = 0;
  Double_t b22 = 0;

  Double_t v1 = h1->GetBinContent(0, coord);
  Double_t v2 = h2->GetBinContent(coord);
  Double_t err1 = h1->GetBinError(coord) * v2;
  Double_t err2 = h2->GetBinError(coord) * v1;
  b22 = v2 * v2;
  err = (err1 * err1 + err2 * err2) / (b22 * b22);

  printf("v1:v2:err1:err2:b22:err = %f:%f:%f:%f:%f:%f\n", v1, v2, err1, err2, b22, err);
  return(TMath::Sqrt(err));
}
