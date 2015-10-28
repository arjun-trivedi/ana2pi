double ye_vgflux(double Ei,double val_W, double val_Q2){
   double mp = 0.938272;
   double a = 0.00729927007299270073;

   double mu = (val_Q2+val_W*val_W-mp*mp)/(2*mp);
   cout <<"mu="<<mu<<endl;
   double Eout = Ei-mu;
   double epilo = 1/(1+2*(val_Q2+mu*mu)/(4*Ei*Eout-val_Q2));
   cout<<"epilo="<<epilo<<endl;
   double t1=(a/(4*TMath::Pi()));
   double t2=(1/(TMath::Power(Ei,2)*TMath::Power(mp,2)));
   double t3=((val_W*(val_W-TMath::Power(mp,2)))/((1-epilo)*val_Q2));
   cout <<"t1="<<t1<<endl;
   cout <<"t2="<<t2<<endl;
   cout <<"t3="<<t3<<endl;
   double flux=t1*t2*t3;
   //double flux = (a/(4*TMath::Pi()))*(1/(TMath::Power(Ei,2)*TMath::Power(mp,2)))*((val_W*(val_W-TMath::Power(mp,2)))/((1-epilo)*val_Q2));
   return flux;
}
