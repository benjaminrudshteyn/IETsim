
/* prect  program  written by John Mitchell, based on method of MP
   Generates k points for a centered 2d rectangular lattice.  Output vectors
   are in fractional coordinates of along the primitive axes of
   reciprocal space, which do not coincide with cartesian axes for this lattice 
   The program assumes that a>b in direct space, so arrange your axes thusly.
   Further, the shape of the BZ and irreducible wedge depend on a and b through
   the variables M and N -- for further information, see the program develop-
   ment documentation 
   Modified for EHMACC input format, 6/23/89 -- J.M. 
   Modified to work with bind 03/19/97 -- Greg Landrum 
   */



#include <stdio.h>
main()
  
{
  FILE  *outfile;
  int w,q,j;
  float a,b;   /* lattice parameters, a>b*/
  double M,N;
  double kx,ky,kz,del,max,min;   /* coord. along cartes. axes in units
				    of 4pi/a and 4pi/b */
  double k1,k2,k3;      /* fractional components along r.l. vectors*/
  int progflag;
  
  kz = 0.0;
  k3 = 0.0;    /* 2d problem */
  outfile = fopen("kout", "w");
  printf("This program assumes that the conventional cell has a>b\n");
  printf("If this is not the case, you must switch coordinates to use\n");
  printf("this program.  Additionally, the primitive lattice vectors\n");
  printf("are assumed to point in the cartesian dirctions [a/2,-b/2], \n");
  printf("and [a/2,b/2] to give the proper BZ.\n");
  printf("Enter your axial length a \n");
  scanf("%f",&a);
  printf("Now enter the axial length b \n");
  scanf("%f",&b);
  if (a<b) {
    printf("You have entered a less than b -- the k points will be\n");
    printf("invalid unless you switch your coordinate axes.\n"); 
  }
  else if (a==b){
    printf("You entered values for a centered square lattice.\n");
    printf("This is not a valid lattice.  Use the program\n");
    printf("square.r with your correct lattice vectors.\n");
  }
  else { 
    printf("Enter Your Value of q---make it even\n");
    scanf("%d",&q);
    del = 1.0/q;
    max = (q-1.0)/(2.0*q);
    min = 0.5 - max;
    M = (a*a + b*b)/(4*a*a);
    N = b*b/(a*a);
    j = 0;
    kx = max;
    while (kx>=0) {
      ky = min;
      while (ky>=0 && ky <= M - N*kx+0.000001) {
	w = 4;
	if (kx == 0 || ky == 0)
	  w = 2;
	if (fabs(ky-M+N*kx)<0.000001)
	  w = 2;
	if (kx==0 && ky==0)
	  w=1;
	/* convert to r.l. vector fractional components */
	k1 = kx - ky;
	k2 = kx + ky;
	fprintf(outfile,"%12.8f",k1);
	fprintf(outfile,"%12.8f",k2);
	fprintf(outfile,"%12.8f",k3);
	fprintf(outfile,"%4.0d\n",w);
	++j;
	ky = ky + del; }
      kx = kx - del; }
    fclose(outfile);
    printf("Number of k-points generated = %d\n",j);
  }
}

