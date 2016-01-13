#ifndef WRPR_E_CORR_SUB
#define WRPR_E_CORR_SUB

//! Fortran subroutines
extern"C" {
void e_corr_sub_(float *thetaeld,float* phield,float* pel,float* torcur,int* secte,float* thetaeldnew,float* newpel);
}

//! C++ interface to Fortrain subroutines
inline void e_corr_sub(float theta,float phi,float p,float torcur,int sector,float& theta_corr,float& p_corr ){
  e_corr_sub_(&theta,&phi,&p,&torcur,&sector,&theta_corr,&p_corr);
  return;
}


#endif //WRPR_E_CORR_SUB
