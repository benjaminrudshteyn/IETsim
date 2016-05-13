
/* simple cubic */
/* written by John Mitchell */
/* conversion from PASCAL to C by John Chung */
/* this program calculates cartesian coordinates of k points in fractions
	of 2pi/a--note for this lattice that cartesian and reciprocal
	axes are coincident...the weight associated with point
	= Order[g(k)]/Order[g(k=0)]  */
/* Modified for EHMACC input file format, 6/23/89 -- J.M. 
   Modified to work with bind 03/19/97 -- Greg Landrum */

#include <stdio.h>

main()
{
  
  FILE *outfile;
  int w,q,j;
  double kx,ky,kz,del,max;
  int progflag;
  
  
  
  j = 0;
  outfile = fopen("kout","w");
  /* HG added the next 15 lines*/
  printf("\n");
  printf("You will enter a seed number next.  Here are some helpful hints: \n");
  printf("\n");
  printf("         SEED            NUMBER OF K POINTS\n\n");
  printf("          2                       1\n");
  printf("          4                       4\n");
  printf("          6                      10\n");
  printf("          8                      20\n");
  printf("         10                      35\n");
  printf("         12                      56\n");
  printf("         14                      84\n");
  printf("\n\n In general you will get:\n");
  printf("(n^3)/48 + (n^2)/8 + (n)/6 points, where\n");
  printf("n is your seed number, which must be even.\n");
  printf("\n");                                       
  printf("Enter Your Value of q -- make it even\n");
  scanf("%d",&q);
  del = 1.0/q;
  max = (q-1.0)/(2.0*q);
  ky = max;
  
  while (ky>=0)
      {
	kx = ky;
	while (kx>=0)
	    {
	      kz = kx;
	      while (kz>=0)
		  {                /*assign weights */
		    w=48;
		    if (kx==0 || ky==0 || kz==0)
		      w = 24;
		    if (kx==ky || kx==kz || ky==kz)
		      w = 24;
		    if (kz==0 && kx==ky)
		      w = 12;
		    if (kx==ky && ky==kz)
		      w = 8;
		    if (kx==0 && kz==0)
		      w = 6;
		    if (kx==0 && ky==0 && kz==0)
		      w = 1;
		    fprintf(outfile,"%12.8f",kx);
		    fprintf(outfile,"%12.8f",ky);
		    fprintf(outfile,"%12.8f",kz);
		    fprintf(outfile,"%4.0d\n",w);
		    kz = kz - del;
		    ++j;	
		  }
	      kx = kx - del;
	    }
	ky = ky - del;
      }
  printf("the number of k points generated is...");
  printf("%2.0d",j);
  printf("\n");
}


