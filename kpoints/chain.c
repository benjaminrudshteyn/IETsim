
/* chain   by John Mitchell, 4/29/88.  Generator of special k-points for a
   one dimensional lattice.  Resultant vectors are fractional coordinates of
   the r.l. vector (2pi/a,0,0) 
   Modified for EHMACC format  6/23/89  -- J.M. 
   Modified to work with bind 03/19/97 -- Greg Landrum 
   */

#include <stdio.h>

main ()
  
{
  
  FILE *outfile;
  double k;   /*fractionl coord. of point in units of 2pi/a*/
  double del;
  int q,w,j;      /* number of k points specified by user and associated 				weight*/
  int progflag;
  
  j = 0;
  outfile = fopen("kout","w");
  printf("Enter how many points you want  \n");
  scanf("%d",&q);
  del = 0.5/q;
  k = del/2.0;
  while (k<=0.5) {
    /*assign weight*/
    w = 2;
    if (k<0.0000001 || abs(k-0.5) < 1e-5)
      w = 1;
    fprintf(outfile,"%12.8f",k);
    fprintf(outfile,"%12.8f",0.0);
    fprintf(outfile,"%12.8f",0.0);
    fprintf(outfile,"%4.0d\n",w);
    k = k + del;
    ++j;
  }
  printf("Number of k points generated = %d\n",j);
  fclose(outfile);
}




