#ifndef WRPR_CUT_SC_PD_E16_H
#define WRPR_CUT_SC_PD_E16_H

//! Fortran subroutines
extern"C" {
void myscint_(int* paddle,bool* bad_sc,int* at_mod);
}

//! C++ interface to Fortrain subroutines
inline bool is_scpd_bad_e16(int paddle,int at_mod){
  bool bad_sc=false;	
  myscint_(&paddle,&bad_sc,&at_mod);
  return bad_sc;
}
#endif //WRPR_CUT_SC_PD_E16_H
