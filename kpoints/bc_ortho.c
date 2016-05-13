/* bcc -- by John Mitchell 6/11/88 -- based on method of MP.  A generator
   of special k points for the b.c. orthorhombic lattice.
   The coordinates are fractions of the recip. lattice vectors, which form
   a face centered orthorhombic lattice.  Note that these k point components
   are not cartesian coordinates, even though cartesian coordinates are used
   in the program.  
   Modified to output to the EHMACC format, 6/23/89 -- J.M.
   Modified to exit when bad lattice vectors are entered 4/21/94 Greg Landrum
   Modified to work with bind 03/19/97 -- Greg Landrum 
   */

#include <stdio.h>
#include <math.h>

/* ---------------- global variables ----------------------------- */


double bovasq, covbsq, covasq;   /*squares of axial ratios for defining
				   irreducible wedge boundaries*/

/* --------------------------------------------------------------- */

int weight (kx,ky,kz)
  double kx,ky,kz;
{
  int w;  /* symmetry weight for a given k point */
  
  
  w = 8; /* default general weight */
  if (kz<0.00001 || ky<0.00001 || kx<0.00001)
    w = 4;
  if (fabs(kx-0.5)<0.000001)
    w = 4; 
  if (fabs(kx-(0.25*(1 + 1.0/covasq) - (1.0/covasq)*kz))<0.000001)
    w = 4;
  if (fabs(ky-(0.25*(1 + 1.0/covbsq) - (1.0/covbsq)*kz))<0.000001)
    w = 4;
  if (fabs(kx+(1.0/bovasq)*ky - 0.25*(1.0+(1.0/bovasq)))<0.000001)
    w = 4;
  if (kz<0.000001 && fabs(kx-0.5)<0.000001)
    w = 2;
  if (fabs(kx-0.5)<0.000001 && ky<0.000001)
    w = 2;
  if (kx<0.000001 && fabs(ky-0.25)<0.00001 && fabs(kz-0.25)<0.00001)
    w = 2;
  if (kz<0.000001 && fabs(ky-0.25)<0.00001 && fabs(kx-0.25)<0.00001)
    w = 2;
  if (ky<0.000001 && fabs(kx-0.25)<0.00001 && fabs(kz-0.25)<0.00001)
    w = 2;
  if (kx<0.000001 && ky<0.000001)
    w = 2;
  if (kx<0.000001 && kz<0.000001)
    w = 2;
  if (kz<0.000001 && ky<0.000001)
    w = 2;
  if (fabs(kx-0.25)<0.00001 && fabs(ky-0.25)<0.00001 && fabs(kz-0.25)<0.0001)
    w = 2;
  if (fabs(kx-0.5)<0.000001 && ky<0.000001 && kz<0.000001)
    w = 1;
  return(w);
}


int irred(kx,ky,kz)
  double kx, ky, kz;
{
  int flag;     /* flag = 1 if the k point in question lies in irred. wedge */
  
  flag = 0;   
  /* region 1 */
  if ((kx>0.0 && kx<0.2500001) && (ky>0.0 && (ky<0.2500001-0.25*bovasq*(4*kx-1))) && (kz>0.0 && (kz<0.2500001-covasq*(0.25-kx))))
    flag = 1;
  
  /*  region 2 */
  if ((kx > 0.0 && kx<0.250001) && (ky>=(0.24999+0.25*bovasq*(4*kx-1)) &&
				    (ky<0.250001-0.25*bovasq*(1-4.0*kx))) && (kz>0.0 && 
									      (kz<0.250001-covbsq*(0.25-ky))))
    flag = 1;
  
  /* region 3 */
  if ((kx>0.2499999 &&  kx<0.5000001) && (ky>0.0 && ky<0.250001+0.25*bovasq*
					  (1-4.0*kx)) && (kz>0.0 && kz<0.250001 + covasq*(0.25-kx)))
    flag = 1;
  return (flag);
}





main ()
  
{
  FILE *outfile;
  int wt,q,j;
  double kx,ky,kz;  /*cart. components, units of 4pi/a, 4pi/b, or 4pi/c*/
  double k1,k2,k3;  /* fractional components along r.l. vectors */
  double max, min, del;   
  double a,b,c;     /* axial lengths of unit cell*/
  int flag,progflag;         /* boolean flag */
  
  
  /* ------------------- main program ---------------------- */
  
  outfile = fopen("kout","w");
  printf("Your direct lattice vectors must point in the cartesian ");
  printf("directions: \n");
  printf("[a/2,b/2,-c/2], [-a/2,b/2,c/2], and [a/2,-b/2,c/2]\n");
  printf("\n");
  printf("In addition, you must have a > b > c. If you do not, please\n");
  printf("use program PERMUTE to transform your axes and coordinates.\n");
  
  
  /* loop until the axes are input with a > b > c */
  flag = 0;
  while (! flag){
    printf("Please enter your value for a, b, and c in this order\n");
    scanf("%lf %lf %lf", &a, &b, &c);
    if (a>b && a>c && b>c) 
      flag = 1;
    else{
      printf("Your axes are incorrect. Use PERMUTE if needed.\n");
      exit(-1);
    }
  }
  printf("Please enter your value of q -- make it even \n");
  scanf("%d",&q);
  j = 0;
  del = 1.0/q;
  max = (q-1.0)/(2.0*q);
  min = 0.5 - max;
  kx = max;
  covasq = pow(c/a,2.0);
  bovasq = pow(b/a,2.0);
  covbsq = pow(c/b,2.0);
  /* loop through the entire irreducible wedge */
  while (kx >= 0.0 ) {
    ky = min;
    while (ky<=max)
	{
	  kz = min;
	  while (kz <= max) {
	    /* test to see if the current k point is in the irred. wedge */
	    if (irred (kx,ky,kz))
	      /* if it is in irred. wedge, proceed to assign weight and write it */
		{
		  /*assign proper weights*/
		  wt = weight(kx,ky,kz);
		  /*convert to fract. coord. along r.l. vectors*/
		  k1 = kx + ky - kz;
		  k2 = ky + kz - kx;
		  k3 = kx + kz - ky;
		  /*write to file kout*/
		  fprintf(outfile,"%12.8f",k1);
		  fprintf(outfile,"%12.8f",k2);
		  fprintf(outfile,"%12.8f",k3);
		  fprintf(outfile,"%4.0d\n",wt);
		  ++j;
		}    /*  end the if in irred. routine */
	    kz = kz + del;
	    
	  }              /* end kz loop*/
	  ky = ky + del;
	}   /* end ky loop*/
    kx = kx - del;
  }  /* end kx loop*/
  printf("Number of k points generated = %d\n",j);
}			




