{
//First extract exp pol. obs.
plotr(1);
plotr(5);
gSystem->Exec("makepdf2 1");
gSystem->Exec("makepdf2 5");

//Now see what it looks like in simulation
plotr(1,1);
plotr(5,1);
gSystem->Exec("makepdf2 1 sim");
gSystem->Exec("makepdf2 5 sim");
}
