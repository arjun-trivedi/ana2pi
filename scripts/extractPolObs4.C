{
//First extract exp pol. obs.
plotphi4(1);
plotphi4(5);
gSystem->Exec("makepdf4 1 exp");
gSystem->Exec("makepdf4 5 exp");

//Now see what it looks like in simulation
plotphi4(1,1);
plotphi4(5,1);
gSystem->Exec("makepdf4 1 sim");
gSystem->Exec("makepdf4 5 sim");
}
