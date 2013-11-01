{
  TString gProcOrder_simYield = "eid:efid:qskim:pid:top";
  TString gProcOrder_expYield = "eid:efid:qskim:mom:pid:top";
  TString gProcOrder_mcYield  = "top";

  //mon Mode i.e. monitor cuts as they are applied
  TString gProcOrder_monMode_simYield = "eidmon:efidmon:qskim:pidmon:top";
  TString gProcOrder_monMode_expYield = "eidmon:efidmon:qskim:mom:pidmon:top";

  //mononly Mode i.e. monitor cuts after final event selection
  TString gProcOrder_monOnlyMode_simYield = "eidmon:efidmon:qskim:pidmon:top:eidmononly:efidmononly:pidmononly";
  TString gProcOrder_monOnlyMode_expYield = "eidmon:efidmon:qskim::mom:pidmon:top:eidmononly:efidmononly:pidmononly";
}
