#include "epTree_Looper.h"
#include <cstdio>
#include <iostream>
using namespace std;

int main(int argc,char *argv[]){
  epTree_Looper* m=new epTree_Looper();
  string binw;
  //! get binw
  if (argc<2){
    printf("Enter theta bin width = binw-1 or binw-05\n");
    getline(cin,binw);
  }else if (argc==2){
    binw=argv[1];
  }else{
    printf("Usage:\n");
    printf("./makeTT [binw]\n");
    printf("Optional argument binw=binw-1 or binw-05\n");
    printf("If binw is not entered, user will be prompted to do so by the program\n");
    return 1;
  }
  //! makeTT as per binw
  if (binw=="binw-1") {
    m->fill_hTTnorm(32,14,46);
  }else if (binw=="binw-05"){
    m->fill_hTTnorm(64,14,46);
  }else{
    printf("%s not understood\n",binw.c_str());
    printf("theta bin width = binw-1 or binw-05 only\n");
  }
  return 1;
  
}

