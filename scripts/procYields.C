void procYields(char* infname, TString topList, 
	            Int_t nQ2bins, Float_t Q2l, Float_t Q2h, 
	            Int_t nWbins, Float_t Wl, Float_t Wh,
	            bool useHel=false){

  printf("Going to process %s and produce observables for:\n ",infname);
  printf("[Q2l,Q2h] [Wl,Wh] = [%4.3f,%4.3f] [%4.3f,%4.3f]\n", Q2l, Q2h, Wl, Wh);
  printf("useHel=%s\n", useHel?"true":"false" );
  ProcYields d(infname, topList, nQ2bins, Q2l, Q2h, nWbins, Wl, Wh,useHel);
  d.Proc();
  printf("Done processing and extracting observables for %s\n ",infname);
  return;
}
