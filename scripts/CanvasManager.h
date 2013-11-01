#ifndef CANVASMANAGER_H
#define CANVASMANAGER_H

#include <TCanvas.h>
#include <TObject.h>

#include <vector>
using namespace std;

class CanvasManager{
  public:
    CanvasManager(Int_t cwidth, Int_t cheight);
    ~CanvasManager();

    static const Int_t laptopw = 1366;
    static const Int_t laptoph = 768;
    
    vector<TObject*> v_olist;
    vector<TCanvas*> v_clist;
    Int_t c_ww;
    Int_t c_wh;
    Int_t c_wtopx;
    Int_t c_wtopx0;
    Int_t c_wtopy;
    Int_t c_wtopy1;
    Int_t c_wtopy2;
  
    Int_t cnum; 
    UInt_t cnum_maxx;
    
    void add(TObject* obj);
    void draw(char* option=NULL);
    void save();
    void close();
};
#endif //CANVASMANAGER_H
