//void setupchain_2pisim(bool usePROOF=true, char* nworkers = "workers=6"){
  //if(usePROOF){
    //TProof::Open(nworkers);
  //}  
{
  TProof::Open("workers=6");

  TFileCollection fc("fileList", "", "h10.2pisim.lst");
  TChain *chain = new TChain("h10");
  chain->AddFileInfoList((TCollection*) fc.GetList());

  //if(usePROOF){
    chain->SetProof();
  //}
  //chain->Draw("mctheta");
}
