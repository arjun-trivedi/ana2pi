{
  TCanvas ce("eDC","eDC",400,400);
  ce->Divide(1,2);
  ce->cd(1);
  h10->Draw("vz>>h_vze(100,-40,-10)","id[0]==11 && id==11");
  ce->cd(2);
  h10->Draw("p>>h_pe(100,0,6)","id[0]==11 && id==11");

  TCanvas cp("pDC","pDC",400,400);
  cp->Divide(1,2);
  cp->cd(1);
  h10->Draw("vz>>h_vzp(100,-40,-10)","id[0]==11 && id==2212");
  cp->cd(2);
  h10->Draw("p>>h_pp(100,0,6)","id[0]==11 && id==2212");
}
