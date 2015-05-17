extern"C" {
bool fidu_e_sub_(float *pel,float *thetael,float *phiel,bool *in_fid_reg);
}

bool fidu_e_sub(float pel,float thetael,float phiel,bool in_fid_reg){
  return fidu_e_sub_(&pel,&thetael,&phiel,&in_fid_reg);
}

/* Returns true if particle is within fiducial volume
  + [05-17-15] Currently implemented only for electrons => 'id' is unused
*/
bool Fiducial_e16(int id, float p, float theta, float phi)
{
  bool in_fid_reg=false;
  return fidu_e_sub(p,theta,phi,in_fid_reg);
}
