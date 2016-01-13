#ifndef WRPR_PR
#define WRPR_PR

//! Fortran subroutines
extern"C" {
void preloss_(float *p_me,float* thet,float* pcorr);
}

//! C++ interface to Fortrain subroutines
inline void preloss(float p,float theta,float& p_corr){
  preloss_(&p,&theta,&p_corr);
  return;
}


#endif //WRPR_PR
