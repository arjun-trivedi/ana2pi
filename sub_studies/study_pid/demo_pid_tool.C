{
 gROOT->ProcessLine(".L pid.cpp++");
 Pid pid("exp"); 

 TH2F* dtVp_pre_cut=new TH2F("pre","pre",100,0,5,250,-5,5);
 TH2F* dtVp_pst_cut=new TH2F("pst","pst",100,0,5,250,-5,5);

 for (int i=0;i<100000;i++){
   float p=gRandom->Uniform(0.1,5);
   float dt=gRandom->Uniform(-5,5);
   dtVp_pre_cut->Fill(p,dt);
   if (pid->is_pip(dt,p)){
     dtVp_pst_cut->Fill(p,dt);
   }
 }
 TCanvas *c =new TCanvas();
 c->Divide(1,2);
 c->cd(1);
 dtVp_pre_cut->Draw("colz");
 c->cd(2);
 dtVp_pst_cut->Draw("colz"); 
}
