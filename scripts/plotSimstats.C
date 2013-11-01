void plotSimstats(Int_t Top, Int_t Varset){
  gROOT->ProcessLine(".L plot_simstats.C+");
  plot_simstats(Top, Varset);
}
