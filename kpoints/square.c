
/* square program by John Mitchell, 4/26/88 
   k point generator for 2d square lattices by method of Monkhorst & Pack 
   q controls number of k points generated 
   note that output are fractional coordinates in units of 2pi/a 
   along the cartesian axes, which coincide with the r.l. vectors in this case
   for further details, see folder containing documentation and listings */
/* Modified for EHMACC input file format, 6/23/89 -- J.M. 
   Modified to work with bind 03/19/97 -- Greg Landrum */

#include <stdio.h>
main()
  
{
  FILE  *outfile;
  int w,q,j;
  double kx,ky,kz,del,max;
  int progflag;
  
  
  kz = 0.0;
  outfile = fopen("kout", "w");
  printf("Enter Your Value of q---make it even\n");
  scanf("%d",&q);
  del = 1.0/q;
  max = (q-1.0)/(2.0*q);
  j = 0;
  kx = max;
  while (kx>=0) {
    
    ky = kx;
    
    while (ky>=0) {
      
      w = 8;
      if (kx==ky)
	w=4;
      if (kx==0 || ky==0)
	w=4;
      if (kx==0 && ky==0)
	w=1;
      fprintf(outfile,"%12.8f",kx);
      fprintf(outfile,"%12.8f",ky);
      fprintf(outfile,"%12.8f",kz);
      fprintf(outfile,"%4.0d\n",w);
      
      ++j;
      ky = ky - del; }
    kx = kx - del; }
  fclose(outfile);
  printf("Number of k-points generated = %d\n", j);
}

