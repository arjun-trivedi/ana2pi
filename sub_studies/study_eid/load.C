{
TTree* t=(TTree*)_file0->Get("eid/monitor/t");

gStyle->SetOptStat("neiou");

TCut pre_cut("etot>0.001 && ec_ei>0.001 && ec_eo>0.001");

TCut s1("sector==1");
TCut s2("sector==2");
TCut s3("sector==3");
TCut s4("sector==4");
TCut s5("sector==5");
TCut s6("sector==6");

TCut noMIP("ec_ei>0.06");
TCut noMIP1("ec_ei>0.1");
TCut noMIP2("ec_ei>0.15");
TCut noMIP3("ec_ei>0.25");

TCut MIP("ec_ei<0.06");
TCut MIP_noH("ec_ei<0.06&&ec_eo<0.1");
TCut MIP1("ec_ei<0.03");
TCut MIP1_noH("ec_ei<0.03&&ec_eo<0.1");

TCut Wlt2("W<2");
TCut Wgt2("W>2");

TString eoVei="ec_eo:ec_ei>>heoVei(150,0,1,150,0,1)";
TString etotVp="etot:p>>hetotVp(160,0,6,150,0,2)";
TString sfVp="etot/p:p>>hsfVp(160,0,6,150,0,0.5)";
TString nphe="nphe>>hnphe(100,0,300)";
TString ccthetaVccsegm="cc_theta:cc_segm>>hccthetaVccsegm(20,0,20,100,0,50)";
TString sfVU="etot/p:ecU>>hsfVU(160,0,450,150,0,0.5)";
TString sfVV="etot/p:ecV>>hsfVV(160,0,450,150,0,0.5)";
TString sfVW="etot/p:ecW>>hsfVW(160,0,450,150,0,0.5)";

TCanvas ceoVei;
t->Draw(eoVei,pre_cut,"colz",1000000);

TCanvas csfVp;
t->Draw(sfVp,pre_cut,"colz",1000000);
}
