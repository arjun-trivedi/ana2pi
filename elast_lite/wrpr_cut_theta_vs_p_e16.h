#ifndef WRPR_CUT_THETA_VS_P_E16_H
#define WRPR_CUT_THETA_VS_P_E16_H

//! Fortran subroutines
extern"C" {
void theta_vs_p_el_(int* sect,float* theta,float* p,bool* good);
}
extern"C" {
void theta_vs_p_pr_(int* sect,float* theta,float* p,bool* good);
}
extern"C" {
void theta_vs_p_pip_(int* sect,float* theta,float* p,bool* good);
}


//! C++ interface to Fortrain subroutines
inline bool theta_vs_p_e16_el(int sect,float theta,float p){
  bool good=true;	
  theta_vs_p_el_(&sect,&theta,&p,&good);
  return good;
}
inline bool theta_vs_p_e16_pr(int sect,float theta,float p){
  bool good=true;	
  theta_vs_p_pr_(&sect,&theta,&p,&good);
  return good;
}
inline bool theta_vs_p_e16_pip(int sect,float theta,float p){
  bool good=true;	
  theta_vs_p_pip_(&sect,&theta,&p,&good);
  return good;
}


#endif //WRPR_CUT_THETA_VS_P_E16_H
