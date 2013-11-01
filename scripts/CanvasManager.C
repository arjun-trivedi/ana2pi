#include "CanvasManager.h"

#include <cstdio>

CanvasManager::CanvasManager(Int_t cwidth, Int_t cheight){
  //printf("CanvasManager constructor\n");
  c_ww = cwidth;
  c_wh = cheight;
  c_wtopx0 = 0;
  c_wtopy1 = 0;
  c_wtopy2 = c_wh;

  cnum = 0;
  cnum_maxx = Int_t(laptopw/cwidth);
  //printf("max = %d\n", cnum_maxx);
}

CanvasManager::~CanvasManager(){
  //printf("CanvasManager destructor\n");
}

void CanvasManager::add(TObject* obj){
  v_olist.push_back(obj);
  return;
}

void CanvasManager::draw(char* option){
  for(UInt_t i = 0; i<v_olist.size(); i++){
        TString opt = option;

        if (v_clist.size() < cnum_maxx){ 
         cnum = v_clist.size();
        }else{
         cnum = v_clist.size() - cnum_maxx;
        }
        //printf("cnum = %d\n", cnum);
        Int_t offsetw = (cnum*c_ww);
        if (v_clist.size() < cnum_maxx){
         c_wtopx = offsetw+c_wtopx0;
         c_wtopy = c_wtopy1; 
        }else{
          c_wtopx = offsetw+c_wtopx0;
          c_wtopy = c_wtopy2;    
        }
        //printf("wtopx, wtopy = %d, %d\n", c_wtopx, c_wtopy);
        v_clist.push_back(new TCanvas(v_olist[i]->GetName(), v_olist[i]->GetTitle(), c_wtopx, c_wtopy,c_ww, c_wh));
        //v_clist.back()->SetLogz(); 
	v_olist[i]->Draw("pads");
       
        if(strcmp(opt,"SetLogz")==0){      
          gPad->Update();
          for (Int_t iPad=1;iPad<30;iPad++){
            if (v_clist.back()->cd(iPad) != NULL) gPad->SetLogz(1);
            else break;
          }
        }
  } 
}

void CanvasManager::save(){
  for(UInt_t i = 0; i<v_clist.size(); i++){
    v_clist[i]->Write();
  }
}

void CanvasManager::close(){
  for(UInt_t i = 0; i<v_clist.size(); i++){
    v_clist[i]->Close();
  }
  v_clist.clear();
}
