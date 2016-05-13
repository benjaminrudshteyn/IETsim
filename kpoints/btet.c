/* btet -- by John Mitchell 4/29/88 -- based on method of MP.  A generator
of special k points for the body centered lattice. The coordinates are fractionsof the reciprocal lattice vectors, which form an fcc lattice (which can, of     course be considered a body centered lattice)
Note that the primitive r.l. vectors do not lie along cartesian axes 
I admit that the programming could be a lot more streamlined -- but this is
a working prototype version.
  Modified for EHMACC format 6/22/89 -- J.M. 
  Modified to work with bind and to fix a seg fault 03/19/97 -- Greg Landrum 
*/

#include <stdio.h>

/* function to return absolute value of argument*/



main ()
{
  FILE *outfile;
  int w,q,j;
  double kx,ky,kz;  /*cartesian components in units of 4pi/a, 4pi/c*/
  double k1,k2,k3;  /* fractional components along r.l. vectors */
  double max, min, del;   
  float a,c;      /* axial lengths */
  double m,n,p,t;   /* constants for a given lattice defining irr. wedge*/
  int progflag;
	  

  outfile = fopen("kout","w");
  printf("Your direct lattice vectors MUST point in the cartesian\n");
  printf("directions: [a/2,a/2,c/2], ");
  printf("[a/2,a/2,-c/2], and [a/2,-a/2,-c/2]\n");
  printf("Please enter your value for length a \n");
  scanf("%f",&a);
  printf("Now enter the length c \n");
  scanf("%f",&c);
  if (a==c) 
    printf("You have a bcc lattice.  Use the appropriate program.\n");
  else {
    printf("Please enter your value of q -- make it even \n");
    scanf("%d",&q);
    j = 0;
    del = 1.0/q;
    max = (q-1.0)/(2.0*q);
    min = 0.5 - max;
    m = (c*c-a*a)/(4*c*c);
    n = (c*c+a*a)/(4*c*c);
    p = (c*c+a*a)/(4*a*a);
    t = (c*c)/(a*a);
    ky = max;
    if (a>c) {
      /* loop through the entire irreducible wedge */
      while (ky >= 0.0 ) {
	kx = min;
	while ((kx<=ky+ 0.000001) && (kx+ky<=0.500001)) {
	  kz = min;
	  while ((p-t*ky)-kz+0.000001 >=0.0)  {
	    /*assign weights*/
	    w = 16; 
	    if (kx<0.000001 || kz<0.00001)
	      w=8;
	    if (fabs(ky-0.25)<1.0e-5 && fabs(kz-.25)<1.0e-5)
	      w = 8;
	    if (fabs(kx-ky)<1e-5)
	      w = 8;
	    if (fabs(kz-0.5)<1e-5)
	      w = 8;
	    if(fabs(kx-ky)<1e-5 && kz<1e-5)
	      w = 4;
	    if(fabs(kz-.25)<1e-5&&fabs(ky-.25)<1e-5&&kx<1e-5)
	      w = 4;
	    if(kz<1e-5 && fabs(kx+ky-.5)<1e-5)
	      w = 4;
	    if(fabs(kx-.25)<1e-5 && fabs(ky-.25)<1e-5)
	      w = 4;
	    if(kx<1e-5 && kz<1e-5)
	      w = 4;
	    if(fabs(ky-.25)<1e-5&&fabs(kx-.25)<1e-5&&kz<1e-5)
	      w = 2;
	    if(kx<1e-5 && fabs(ky-.5)<1e-5)
	      w = 2;
	    if(fabs(kx-ky)<1e-5&&fabs(kx-kz)<1e-5&&fabs(kx-0.25)<1e-5)
	      w = 2;
	    if(kx<1e-5 && ky<1e-5)
	      w = 2;
	    if(fabs(ky-.5)<1e-5&&kx<1e-5&&kz<1e-5)
	      w = 2;  
	    /*convert to fract. coord. along r.l. vectors*/
	    k1 = kx + ky + kz;
	    k2 = kx + ky - kz;
	    k3 = kx - ky - kz;
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
    }
    else   /* now c>a */ {
      while (ky > n+1e-5) {
	ky = ky - del;
      }
      while (ky > 0.0) {
	kx = min;
	while (kx<ky+1e-5 && kx+ky<0.5000001) {
	  kz = min;
	  while ((kz<0.500001 && ky<m+1e-5)||(ky>m-1e-5
					      && p-t*ky-kz+1e-5 >=0.0)) {
	    /*assign weights*/
	    w = 16;
	    if(kx<1e-5||kz<1e-5)
	      w = 8;
	    if(fabs(ky-kx)<1e-5||fabs(kz-.5)<1e-5)
	      w = 8;
	    if(kx<1e-5 && kz<1e-5)
	      w = 4;
	    if(fabs(kx-ky)<1e-5 && (fabs(kz-.5)<1e-5 
				    || kz<1e-5))
	      w = 4;
	    if(kx<1e-5 && fabs(kz-0.5)<1e-5)
	      w = 4;
	    if(kx<1e-5&&fabs(kz-.25)<1e-5&&fabs(ky-				   	           0.25)<1e-5)
	      w = 4;
	    if(fabs(kx-ky)<1e-5&&fabs(kz-p-t*ky)<1e-5)						w = 4;
	    if(fabs(kx+ky-.5)<1e-5)
	      w = 4;
	    if(fabs(kx+ky-.5)<1e-5&&kz<1e-5)
	      w = 2;
	    if(fabs(kx-.25)<1e-5&&fabs(ky-.25)<1e-5 &&
	       fabs(kz-0.25)<1e-5)
	      w = 2;
	    /*convert to r.l. coords*/
	    k1 = kx + ky + kz;
	    k2 = kx + ky - kz;
	    k3 = kx - ky - kz;
	    fprintf(outfile,"%12.8f",k1);
	    fprintf(outfile,"%12.8f",k2);
	    fprintf(outfile,"%12.8f",k3);
	    fprintf(outfile,"%4.0d\n",w);
	    ++j;
	    kz = kz + del;
	  }  /*end kz loop*/
	  kx = kx + del;
	}    /*end kx loop*/
	ky = ky - del;
      }     /*end ky loop*/
    }     /* end else c>a portion */					
    printf("Number of k points generated = %d\n",j);
  }
}			


	

