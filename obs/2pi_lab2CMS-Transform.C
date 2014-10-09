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
}

