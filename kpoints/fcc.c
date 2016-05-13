/* fcc -- by John Mitchell 4/26/88 -- based on method of MP.  A generator
   of special k points for the fcc lattice.  The coordinates are fractions
   of the reciprocal lattice vectors, which form a bcc lattice.  Note that
   these are not cartesian coordinates ! 
   Modified for EHMACC input file format  6/23/89 -- J.M. 
   Modified to work with bind 03/19/97 -- Greg Landrum 
   */

#include <stdio.h>



main ()
{
  FILE *outfile;
  int w,q,j;
  double kx,ky,kz;  /*cartesian components in units of 4pi/a*/
  double k1,k2,k3;  /* fractional components along r.l. vectors */
  double max, min, del;   
  int progflag;
  
  outfile = fopen("kout","w");
  printf("Your direct lattice vectors must point in the cartesian ");
  printf("directions: [a/2,0,a/2],\n");
  printf("[a/2,a/2,0], and [0,a/2,a/2]\n");
  printf("Please enter your value of q -- make it a multiple of 4\n");
  scanf("%d",&q);
  j = 0;
  del = 1.0/q;
  max = (q-1.0)/(2.0*q);
  min = 0.5 - max;
  ky = max;
  /* loop through the entire irreducible wedge */
  while (ky >= 0.0 ) {
    kx = min;
    while (kx >=0 && (kx<=ky+ 0.000001) && (kx + ky <= 0.7500001)) {
      kz = min;
      while (kz >= 0 && (kz <= kx+0.00001) && (kx+ky+kz <=    			 0.7500001)) {
	/*assign proper weights*/
	w = 48; /*default*/
	if (fabs(kx-ky)<0.000001 || fabs(ky-0.5)<0.00001)
	  w = 24;
	if (kz < 0.0000001) 
	  w=24;
	if (fabs(kx+ky+kz-0.75)<0.000001 || 
	    fabs(kx-kz)<0.0001)
	  w = 24;
	if (fabs(kx-0.25)<0.000001 &&  
	    fabs(kz+ky-0.5)<0.00001)
	  w=24;
	if (fabs(ky-0.5)<0.00001 && kz<0.00001) 
	  w=12;
	if (fabs(ky-0.5)<0.00001 && fabs(kx-kz)<0.00001)
	  w=12;
	if (kz<0.00001 && fabs(kx-0.375)<0.000001 && 
	    fabs(ky-0.375)<0.000001)
	  w=12;
	if (fabs(kx-ky)<0.00001 && kz<0.00001)
	  w=12;
	if (fabs(kx-ky)<0.00001 && fabs(ky-kz)<0.00001)
	  w=8;
	if (kx<0.000001 && kz<0.000001)
	  w=6;
	if (kz<0.000001 && fabs(kx-0.25)<0.000001 
	    && fabs(ky-0.5)<0.00001)
	  w=6;
	if (fabs(kx-0.25)<0.00001 && fabs(ky-0.25)<0.0001					   && fabs(kz-0.25)<0.00001)
	  w=4;
	if (kx<0.00001 &&  kz<0.000001 &&  fabs(ky-0.5)< 				   0.000001)
	  w=3;
	/*convert to fract. coord. along r.l. vectors*/
	k1 = kx + ky;
	k2 = ky + kz;
	k3 = kx + kz;
	/*write to file kout*/
	fprintf(outfile,"%12.8f",k1);
	fprintf(outfile,"%12.8f",k2);
	fprintf(outfile,"%12.8f",k3);
	fprintf(outfile,"%4.0d\n",w);
	kz = kz + del;
	++j;
      }              /* end kz loop*/
      kx = kx + del;
    }   /* end kx loop*/
    ky = ky - del;
  }  /* end ky loop*/
  printf("Number of k points generated = %d\n",j);
}			


	

