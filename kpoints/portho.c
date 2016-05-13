
/* Written by John Mitchell 04/27/88
   converted to C from PASCAL by John Chung
   
   this program calculates the cartesian coordinates of k points
   in fractions of 2pi/a, 2pi/b, 2pi/c -- for this lattice
   the cartesian and reciprocal axes are coincident
   
   weight associated with point = Order[g(k)]/order[g(k=0)]
   
   ******************************************************** */
/* Modified for EHMACC input file format, 6/23/89 -- J.M. 
   Modified to work with bind 03/19/97 -- Greg Landrum */

#include <stdio.h>

main()
  
{
  
  FILE *outfile;
  int w,q,j;
  double kx,ky,kz,del,max;
  int progflag;
  
  outfile = fopen("kout","w");
  /* HG added the next 15 lines*/
  printf("\n\n\n");
  printf("You will enter a seed number next.  Here are some helpful hints: \n");
  printf("\n");
  printf("         SEED            NUMBER OF K POINTS\n\n");
  printf("          2                       1\n");
  printf("          4                       8\n");
  printf("          6                      27\n");
  printf("          8                      64\n");
  printf("         10                     125\n");
  printf("         12                     216\n");
  printf("         14                     343\n");
  printf("\n\n In general you will get:\n");
  printf("(n^3)/8 points, where\n");
  printf("n is your seed number, which must be even.\n");
  printf("\n");
  printf("Enter your value of q--make it even!\n");
  scanf("%d",&q);
  j = 0;
  del = 1.0/q;
  max = (q-1.0)/(2.0*q);
  kx = max;
  while (kx>=0) {
    ky = max;
    while (ky>=0) {
      kz = max;
      while (kz>=0) {
	w = 8;
	if (kx==0 || ky==0 || kz==0)
	  w=4;
	if (ky==0 && kz==0)
	  w=4;
	if (kx==0 && ky==0 && kz==0)
	  w=1;
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
  printf("%2.0d",j);
  printf("k points were generated!\n");
}

