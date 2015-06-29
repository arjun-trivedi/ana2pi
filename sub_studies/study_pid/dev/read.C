void read(){
 FILE *fp = fopen("/home/trivedia/CLAS/workspace/ana2pi/sub_studies/study_pid/dtVp_cuts/Exp/p/cutpars.txt", "r");
 //FILE *fp = fopen("test.txt", "r");
 if (fp) {
   fprintf(stdout, "Reading file\n");
   char str[200];
   char n[20];
   float a,b,c,d;
   //int a,b,c,d;
   while (fgets(str, sizeof(str), fp)) {
     cout<<str<<endl;
     //sscanf(str, "%s %d %d %d %d", n, &a, &b,&c,&d);
     //fprintf(stdout, "%s %d %d %d %d\n", n,a,b,c,d);
     sscanf(str, "%s %f %f %f %f", n,&a, &b,&c,&d);
     fprintf(stdout, "%s %f %f %f %f\n",n,a,b,c,d);
   }
 }else{
   printf("File cannot be read!\n");
 }
 fclose(fp);
 return;
}

