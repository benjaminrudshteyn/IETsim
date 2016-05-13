
/* simple tetragonal */
/* written by John Mitchell */
/* this program calculates cartesian coordinates of k points in fractions
   of 2pi/a, 2pi/c--note for this lattice that cartesian and reciprocal
   axes are coincident...the weight associated with point
   = Order[g(k)]/Order[g(k=0)]  */
/* Modified for EHMACC input file format, 6/23/89 -- J.M. 
   Modified to work with bind 03/19/97 -- Greg Landrum */

#include <stdio.h>

main()
{
  
  FILE *outfile;
  int w,q,j=0;
  double kx,ky,kz,del,max;
  int progflag;
  
  outfile = fopen("kout","w");
  printf("Enter Your Value of q -- make it even\n");
  scanf("%d",&q);
  del = 1.0/q;
  max = (q-1.0)/(2.0*q);
  kx = max;
  
  while (kx>=0)
      {
	ky = kx;
	while (ky>=0)
	    {
	      kz = max;
	      while (kz>=0)
		  {                /*assign weighta */
		    w=16;
		    if (kx==0 || ky==0 || kz==0)
		      w = 8;
		    if (kx==ky)
		      w = 8;
		    if (kz==0 && kx==ky)
		      w = 4;
		    if (kz==0  && ky==0)
		      w = 4;
		    if (kx==0 && kz==0 && ky==0)
		      w = 1;
		    fprintf(outfile,"%12.8f",kx);
		    fprintf(outfile,"%12.8f",ky);
		    fprintf(outfile,"%12.8f",kz);
		    fprintf(outfile,"%4.0d\n",w);
		    kz = kz - del;
		    ++j;	
		  }
	      ky = ky - del;
	    }
	kx = kx - del;
      }
  printf("the number of k points generated is...");
  printf("%2.0d",j);
  printf("\n");
}

	
