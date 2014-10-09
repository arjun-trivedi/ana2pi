// [10-09-14] Here, I began to note how the CMS system 
// defined by Evan(ep) and Gleb(gf) are different
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TRotation.h>
int run(){
  //! Set up initial states
  float MASS_P = 0.93827203;
  float MASS_E = 0.000511;
  float E1F_P=5.499;
  TLorentzVector lvP0(0,0,0,MASS_P);
  TLorentzVector lvE0(0,0,E1F_P,TMath::Sqrt(E1F_P*E1F_P+MASS_E*MASS_E));

  //! "Simulate" E1
  TLorentzVector lvE1(0,0,0,0);
  TVector3 vE1(0.62,0.62,3);
  lvE1.SetVectM(vE1,0.000511);

  //! Obtain CMS system basis vectors uz and ux 
  //! + uz is set up direction of virtual photon
  //! + ux is set by E0 and E1
  //! 	+ ep and gf differ in their calculated ux
  TLorentzVector lvQ(0,0,0,0);
  TLorentzVector lvW(0,0,0,0);
  lvQ=lvE0-lvE1;
  lvW=lvQ+lvP0;
  TVector3 uz=lvQ.Vect().Unit();
  TVector3 ux_ep=(lvE0.Vect().Cross(lvE1.Vect())).Unit();
  ux_ep.Rotate(-TMath::Pi()/2,uz);
  TVector3 ux_gf=(lvE1.Vect().Cross(lvE0.Vect())).Unit();
  ux_gf.Rotate(3*TMath::Pi()/2,uz);

  
  //! Calculate Rotation(r3) to transform lab frame vectors to CMS
  TRotation r3_ep;
  r3_ep.SetZAxis(uz,ux_ep).Invert();
  TRotation r3_gf;
  r3_gf.SetZAxis(uz,ux_gf).Invert();

  //! Calculate Evan's Transformation which is a LorentzRotation
  //! and therefore done as a single operation
  //! +  Transformation = LorentzRotation: r4=r3*boost
  TVector3 boost(-1*lvW.BoostVector());
  TLorentzRotation r4(r3_ep); //*_boost);
  r4*=boost;

  //! Calculate Gleb's Transformation which is a two step process, that of 
  //! a spatial rotation followed by a boost:
  //! 	 1. r3
  //!	 2. boost
  float Q2=-lvQ.Mag2();
  float beta=TMath::Sqrt(lvQ[3]*lvQ[3]+Q2)/(lvQ[3]+MASS_P);

  //! Perform Evan's Transformation on lvQ and lvE1
  TLorentzVector lvQ_CMS_ep(0,0,0,0);
  TLorentzVector lvE1_CMS_ep(0,0,0,0);
  lvQ_CMS_ep=lvQ;
  lvE1_CMS_ep=lvE1;
  lvQ_CMS_ep.Transform(r4);
  lvE1_CMS_ep.Transform(r4);
  
  //! Perform Gleb's Transformation on lvQ and lvE1
  TLorentzVector lvQ_CMS_gf(0,0,0,0);
  TLorentzVector lvE1_CMS_gf(0,0,0,0);
  lvQ_CMS_gf=lvQ;
  lvE1_CMS_gf=lvE1;
  lvQ_CMS_gf.Transform(r3_gf);
  lvQ_CMS_gf.Boost(0,0,-beta);
  lvE1_CMS_gf.Transform(r3_gf);
  lvE1_CMS_gf.Boost(0,0,-beta);
  

  printf("uz=(%f,%f,%f)\n",uz.X(),uz.Y(),uz.Z());
  printf("ux_ep=(%f,%f,%f)\n",ux_ep.X(),ux_ep.Y(),ux_ep.Z());
  printf("ux_gf=(%f,%f,%f)\n",ux_gf.X(),ux_gf.Y(),ux_gf.Z());

  printf("r3_ep:\n");
  printf("%f,%f,%f\n",r3_ep.XX(),r3_ep.XY(),r3_ep.XZ());
  printf("%f,%f,%f\n",r3_ep.YX(),r3_ep.YY(),r3_ep.YZ());
  printf("%f,%f,%f\n",r3_ep.ZX(),r3_ep.ZY(),r3_ep.ZZ());
  printf("r3_gf:\n");
  printf("%f,%f,%f\n",r3_gf.XX(),r3_gf.XY(),r3_gf.XZ());
  printf("%f,%f,%f\n",r3_gf.YX(),r3_gf.YY(),r3_gf.YZ());
  printf("%f,%f,%f\n",r3_gf.ZX(),r3_gf.ZY(),r3_gf.ZZ());

  printf("ep boost vector magnitude=%f\n",boost.Mag());
  printf("gleb boost's beta_z=%f\n",beta);

  printf("ep lvQ_CMS=(%f,%f,%f,%f)\n",lvQ_CMS_ep.X(),lvQ_CMS_ep.Y(),lvQ_CMS_ep.Z(),lvQ_CMS_ep.T());
  printf("gf lvQ_CMS=(%f,%f,%f,%f)\n",lvQ_CMS_gf.X(),lvQ_CMS_gf.Y(),lvQ_CMS_gf.Z(),lvQ_CMS_gf.T());
  printf("ep lvE1_CMS=(%f,%f,%f,%f)\n",lvE1_CMS_ep.X(),lvE1_CMS_ep.Y(),lvE1_CMS_ep.Z(),lvE1_CMS_ep.T());
  printf("gf lvE1_CMS=(%f,%f,%f,%f)\n",lvE1_CMS_gf.X(),lvE1_CMS_gf.Y(),lvE1_CMS_gf.Z(),lvE1_CMS_gf.T());

  printf("dbg=%f,%f,%f\n",lvQ_CMS_ep.Mag(), lvQ_CMS_gf.Mag(), -lvQ.Mag2());
  printf("dbg=%f,%f,%f\n",lvE1_CMS_ep.Mag(),lvE1_CMS_gf.Mag(),lvE1.Mag());
  printf("W=%f",lvW.Mag());
  
  return 0;
}

