#ifndef WRPR_VERTEX_E16
#define WRPR_VERTEX_E16

//! Fortran subroutines
extern"C" {
void vertex_e16_(float *px,float* py,float* pz,float* vx0,float* vy0,float* vz0,float* vx,float* vy,float* vz);
}

//! C++ interface to Fortrain subroutines
inline void vertex_e16(float px,float py,float pz,float vx0,float vy0,float vz0,float& vx,float& vy,float& vz){
  vertex_e16_(&px,&py,&pz,&vx0,&vy0,&vz0,&vx,&vy,&vz);
  return;
}

#endif //WRPR_VERTEX_E16
