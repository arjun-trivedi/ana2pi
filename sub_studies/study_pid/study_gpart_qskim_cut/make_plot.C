{
TFile* f=TFile::Open("$D2PIDIR_EXP/data_pid_gpart_qskim_cut/dpid_test.root");
TTree* t=(TTree*)f->Get("pid/cut/t");
gROOT->ProcessLine(".L looper.C+");
looper m(t);
m.Loop();
}
