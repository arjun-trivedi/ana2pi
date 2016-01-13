#ifndef WRPR_PI_CORR_SUB
#define WRPR_PI_CORR_SUB

//! Fortran subroutines
extern"C" {
void pi_corr_sub_(float *thetapid,float* phipid,float* mom_pip,float* torcur,int* secth,float* nthetapidnew,float* newpip);
}

//! C++ interface to Fortrain subroutines
inline void pi_corr_sub(float theta,float phi,float p,float torcur,int sector,float& theta_corr,float& p_corr ){
  pi_corr_sub_(&theta,&phi,&p,&torcur,&sector,&theta_corr,&p_corr);
  return;
}


#endif //WRPR_PI_CORR_SUB
