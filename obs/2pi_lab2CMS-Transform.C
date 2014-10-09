// [10-09-14] Here, I began to note how the CMS system 
// defined by Evan(ep) and Gleb(gf) are different

{
  //! Set up initial states
  float MASS_P = 0.93827203;
  float MASS_E = 0.000511;
  float E1F_P=5.499;
  TLorentzVector lvP0(0,0,0,MASS_P);
  TLorentzVector lvE0(0,0,E1F_P,sqrt(E1F_P**2+MASS_E**2));

  //! "Simulate" E1
  TLorentzVector lvE1(0,0,0,0);
  TVector3 vE1(1.62,1.62,5);
  lvE1.SetVectM(vE1,0.000511);

  //! Obtain CMS system basis vectors uz and ux 
  //! + uz is set up direction of virtual photon
  //! + ux is set by E0 and E1
  //! 	+ ep and gf differ in their calculated ux
  lvQ=lvE0-lvE1;
  lvW=lvQ+lvP0;
  uz=lvQ.Vect().Unit();
  ux_ep=(lvE0.Vect().Cross(lvE1.Vect())).Unit();
  ux_ep.Rotate(-TMath::Pi()/2,uz);
  ux_gf=(lvE1.Vect().Cross(lvE0.Vect())).Unit();
  ux_gf.Rotate(3*TMath::Pi()/2,uz);

  r3_ep=TRotation();
  r3_ep.SetZAxis(uz,ux_ep).Invert();
  r3_gf=TRotation();
  r3_gf.SetZAxis(uz,ux_gf).Invert();

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
}

