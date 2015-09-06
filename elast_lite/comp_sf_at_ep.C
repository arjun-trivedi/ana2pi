{
//! sector 1 only
TF1* l_at_exp=new TF1("at","pol3",0,10);
TF1* h_at_exp=new TF1("at","pol3",0,10);
l_at_exp->SetParameters(0.11831,0.08646,-0.01782,0.00148);
h_at_exp->SetParameters(0.40040,-0.02844,0.01311,-0.00150);

TF1* l_at_sim=new TF1("at","pol3",0,10);
TF1* h_at_sim=new TF1("at","pol3",0,10);
l_at_sim->SetParameters(0.11444,0.05785,-0.01546,0.00152);
h_at_sim->SetParameters(0.32710,-0.02836,0.00567,-0.00049);

TF1* l_ep_exp=new TF1("ep","pol3",0,10);
TF1* h_ep_exp=new TF1("ep","pol3",0,10);
l_ep_exp->SetParameters(0.1,0.12,-0.035,0.0038);
h_ep_exp->SetParameters(0.4,-0.037,0.019,-0.0023);
l_ep_exp->SetLineColor(kBlack);
h_ep_exp->SetLineColor(kBlack);

TF1* l_ep_sim=new TF1("ep","pol3",0,10);
TF1* h_ep_sim=new TF1("ep","pol3",0,10);
l_ep_sim->SetParameters(0.13,0.043,-0.008,0.00021);
h_ep_sim->SetParameters(0.31,-0.013,-0.0016,0.00072);
l_ep_sim->SetLineColor(kBlack);
h_ep_sim->SetLineColor(kBlack);


TCanvas* c=new TCanvas();
c->Divide(1,2);
c->cd(1);
TH1F* h=new TH1F("exp","exp",100,0,10);
h->SetMaximum(0.5);
h->Draw();
l_at_exp->Draw("same");
h_at_exp->Draw("same");
l_ep_exp->Draw("same");
h_ep_exp->Draw("same");
c->cd(2);
TH1F* h=new TH1F("sim","sim",100,0,10);
h->SetMaximum(0.5);
h->Draw();
l_at_sim->Draw("same");
h_at_sim->Draw("same");
l_ep_sim->Draw("same");
h_ep_sim->Draw("same");
}
