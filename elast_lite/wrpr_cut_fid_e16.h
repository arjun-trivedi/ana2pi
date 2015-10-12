#ifndef WRPR_CUT_FID_E16_H
#define WRPR_CUT_FID_E16_H

extern"C" {
bool fidu_e_sub_(float *pel,float *thetael,float *phiel,bool *in_fid_reg);
}
extern"C" {
bool hadronfid_(bool *result,float *theta,float *phi,int *s);
}

bool fidu_e_sub(float pel,float thetael,float phiel,bool in_fid_reg){
  return fidu_e_sub_(&pel,&thetael,&phiel,&in_fid_reg);
}
bool hadronfid(bool result,float theta,float phi,int s){
  return hadronfid_(&result,&theta,&phi,&s);
}

/* Returns true if electron is within fiducial volume
  + 'id' is unused
*/
bool Fiducial_e16_elctrn(int id, float p, float theta, float phi)
{
  bool in_fid_reg=false;
  return fidu_e_sub(p,theta,phi,in_fid_reg);
}

bool Fiducial_e16_hdrn(float theta, float phi,int sector)
{
  bool in_fid_reg=false;
  return hadronfid(in_fid_reg,theta,phi,sector);
}

#endif //WRPR_CUT_FID_E16_H
