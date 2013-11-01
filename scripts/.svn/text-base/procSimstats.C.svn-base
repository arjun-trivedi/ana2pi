void procSimstats(TString sim){
  cout << ".L myTHnTool.C+ ..." << endl;
  gROOT->ProcessLine(".L myTHnTool.C+");
  gInterpreter->AddIncludePath("/home/trivedia/ongoing/macros/");
  cout << ".L proc_simstats.C+ ..." << endl;
  gROOT->ProcessLine(".L proc_simstats.C+");
  //proc_simstats(sim,"at");
  proc_simstats(sim,"vm");
  proc_simstats(sim,"vm",POS);
  proc_simstats(sim,"vm",NEG);
  cout <<"procSimstats::Done"<<endl;
  return;
}
