void runSelector(Long64_t nentries=1000000000, TString fin="", TString fout = "", TString procorder = "") {
 
 /* *** determine fin and setup TChain *** */
 TChain* chain;

 if (fin.Contains(".root")) {
   chain = new TChain("h10");
   chain->Add(fin);
 }else if(fin.Contains(".lst")){
   if(fin.Contains("h10skim")){
     chain = new TChain("skim/h10");
     TFileCollection fc("fileList", "", fin);
     chain->AddFileInfoList((TCollection*) fc.GetList());
   }else{
     chain = new TChain("h10");
     TFileCollection fc("fileList", "", fin);
     chain->AddFileInfoList((TCollection*) fc.GetList());
   }   
 }else{
   Info("runSelector", "Could not determine fin!");
   return;
 }
 //Info("runSelector", "Total number of entries in chain = %d", chain->GetEntries());
  
 /* *** determine H10TYPE = (exp:dtype:skim) *** */
 TString h10type;
 TString h10exp, h10dtype, h10skim;
 Bool_t is_h10e1f, is_h10e16, is_h10exp, is_h10sim = kFALSE;
  
 if      (chain->GetBranch("evthel"))  h10type = "e1f:exp";
 else if (chain->GetBranch("nschit"))  h10type = "e16:exp";
 else if (chain->GetBranch("mcnentr")) h10type = "e1f:sim";
 if (fin.Contains("h10skim")) h10type += ":skim";
  
 TObjArray *h10type_tokens = h10type.Tokenize(":");
 if(h10type_tokens->GetEntries() >= 2){	
   h10exp   = h10type_tokens->At(0)->GetName();
   h10dtype = h10type_tokens->At(1)->GetName();
   if (h10type_tokens->GetEntries() == 3) h10skim = h10type_tokens->At(2)->GetName();
	
   if(h10exp.EqualTo("e1f")) is_h10e1f = kTRUE;
   else if (h10exp.EqualTo("e16")) is_h10e16 = kTRUE;
   else {
     Info("runSelector", "Could not determine h10type.experiment!\n");
     return;
   }
	
   if (h10dtype.EqualTo("exp"))      is_h10exp = kTRUE;
   else if (h10dtype.EqualTo("sim")) is_h10sim = kTRUE;
   else {
     Info("runSelector", "Could not determine h10type.dtype!\n");
     return;
   }
   Info("runSelector","Going to set up SelH10 with h10type: %s:%s:%s\n", h10exp.Data(), h10dtype.Data(), h10skim.Data());
 }else{
  Info("runSelector", "h10type identifier does not have enough information! Tokens < 2\n");
  return;
 }
 
 /* *** setup PROOF *** */
 bool usePROOF=false;
 TString host = gSystem->Getenv("HOST");
 if ( !fin.Contains(".root") && host.EqualTo("gothe14")) usePROOF=true; 
 if (usePROOF){
   TProof::Open("workers=4");
   gProof->Load("/home/trivedia/CLAS/workspace/ana2pi/libAna2pi.so");
   gProof->SetParameter("PROOF_CacheSize", 6000000000);
   chain->SetProof();
 }else{
   gSystem->Load("/home/trivedia/CLAS/workspace/ana2pi/libAna2pi.so");
 }
 
 /* *** set up SelH10 *** */
 SelH10 *s = new SelH10();
 if(usePROOF){
   TList* inputList = new TList();
   inputList->Add(new TNamed("PROOF_OUTPUTFILE", fout.Data()));
   inputList->Add(new TNamed("H10TYPE", h10type.Data()));
   s->SetInputList(inputList);
 }else{
   TList* inputList = new TList();
   inputList->Add(new TNamed("OUTPUTFILE", fout.Data()));
   inputList->Add(new TNamed("H10TYPE", h10type.Data()));
   s->SetInputList(inputList);
 }
 
 Info("runSelector", "Going to run SelH10:\n");
 Info("runSelector", "fin =  %s", fin.Data());
 Info("runSelector", "fout =  %s", fout.Data());
 Info("runSelector", "procorder =  %s", procorder.Data());
 Info("runSelector", "h10type =  %s", h10type.Data()); 
 
 chain->Process(s,procorder, nentries); 
  
 Info("runSelector", "SUCCESS!\n");
}

