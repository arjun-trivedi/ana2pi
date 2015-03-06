extern "C"{
  float elasrad_( float* es, float* theta_d, float* t, float* wcut );
  float elas_( float* eb, float* theta);
}

float diffxsec_elasrad(float es, float theta, float t, float wcut){
  return elasrad_(&es, &theta, &t, &wcut );
}

float diffxsec_elas(float eb, float theta){
  return elas_(&eb, &theta);
}