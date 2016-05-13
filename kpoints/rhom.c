/* rhom.c by John Mitchell 5/11/88 -- program to calculate special k points
   for a rhombohedral lattice -- excluding special cases b.c.c., f.c.c., s.c.
   Uses method of MP.   See documentation folder for derivation of expressions
   for the irreducible wedges of the two BZ's. */
/* Modified for EHMACC input file format, 6/23/89 -- J.M. 
   Modified to work with bind 03/19/97 -- Greg Landrum */

#include <stdio.h>
#include <math.h>


main ()
{
  FILE *outfile;
  double h1,h2,h3;   /* coord. of k points in terms of hexagonal basis */
  double r1,r2,r3;   /* coord. of k points in terms of rhombohedral basis */
  double alphadeg,alpha;      /* angle between basis vectors in direct space */
  double del, max, min,x,y,m;
  int q,w,flag,num;
  int progflag;
  
  outfile = fopen("kout","w");
  printf("Please enter the value of alpha for your direct lattice\n");
  scanf("%lf",&alphadeg);
  alpha = (3.1415926/180.0)*alphadeg;
  while (alphadeg < 0.0 || alphadeg > 120.0) {
    printf("Please enter an angle between 0.0 and 120.0 degrees\n");
    scanf("%lf",&alphadeg);
  }
  if (alphadeg == 0.0) 
    printf("This is a 1-d lattice -- please use appropriate program\n");
  else if (alphadeg == 60.0) 
    printf("This is an f.c.c. lattice -- please use appropriate program\n"); 
  else if (alphadeg == 90.0) 
    printf("This is an s.c. lattice -- please use appropriate program\n");
  else	if (alphadeg == 109.4712)
    printf("This is a b.c.c. lattice -- please use appropriate program\n");
  else if (alphadeg == 120.0)
    printf("This is a 2-d hex lattice -- please use appropriate program\n");
  else {                     /* begin calculation phase */
    printf("Please enter a value for q -- make it even\n");
    scanf("%d",&q);
    del = 1.0/q;
    max = (q-1.0)/(2.0*q);
    min = 0.5 - max;
    x = (3.0 * (1 + 2*cos(alpha)))/(2.0 * (1 - cos(alpha)));  /* sqr(c/a) */
    y = 1.0 / x;
    m = 3 - 4.0*x;
    r1 = max;
    num = 0;
    if (cos(alpha) < 0.0) {
      while (r1 >= -max) {
	r2 = -max;
	while (r2 <= max) {
	  r3 = -max;
	  while (r3 <= max) {
	    h1 = (1.0/3.0) * (2.0 * r1 - r2 - r3);
	    h2 = (1.0/3.0) * (r1 + r2 - 2.0 * r3);
	    h3 = (1.0/3.0) * (r1 + r2 + r3);
	    if (((h1>=0.0 && h1<=1.0/3.0) && (h2>=0.5*h1 && h2<=2.0*h1)||
		 (h1>=1.0/3.0 && h1<=2.0/3.0) &&(h2>=0.5*h1 && h2<=1.0-h1))
		&& (h3>=(m/6.0)*(h1-h2) && h3<=(1.0/6.0 - (2.0/9.0)*x*                          (3*h1-1.0))))  {
	      w = 12;
	      if (h2 == 0.5*h1 || h2 == 2.0*h1) 
		w = 6;
	      if (h1 == h2 && h3 == 0.0) 
		w = 6;
	      if (h1 == 1.0/3.0 && h3 == 1.0/6.0) 
		w = 6;
	      if (h1 == 0.5 && h2 == 0.5 && h3 == 0) 
		w = 3;
	      if (h1 == 1.0/3.0 && h2 == 1.0/6.0 && h3 == 1.0/6.0) 
		w = 3;
	      if (h1 == 0.0 && h2 == 0.0) 
		w = 2;
	      if (h1 == 2.0/3.0 && h2 == 1.0/3.0) 
		w = 2;
	      fprintf(outfile,"%12.8f",r1);
	      fprintf(outfile,"%12.8f",r2);
	      fprintf(outfile,"%12.8f",r3);
	      fprintf(outfile,"%4.0d\n",w);
	      num = num + 1;
	    }
	    r3 = r3 + del;
	  }
	  r2 = r2 + del;
	}
	r1 = r1 - del;
      }
    }
    else   {       /* alpha < 90.0 deg. */
      while (r1>=-max) {            
	r2 = -max;
	while (r2<=max) {
	  r3 = -max;
	  while (r3 <= max) {
	    h1 = (1.0/3.0)*(2.0*r1 - r2 - r3);
	    h2 = (1.0/3.0)*(r1 + r2 - 2.0*r3);
	    h3 = (1.0/3.0)*(r1 + r2 + r3);
	    flag = 0;
	    /*region A*/   if ((((h1>=0.0&&h1<=(1.0/6.0 - 0.25*y)) && (h2>=0.5*h1 &&
								       h2<=2.0*h1)) || ((h1>=(1.0/6.0-0.25*y) && h1<=(1.0/3.0 -
														      0.5*y)) && (h2>=0.5*h1 && h2<=(1.0/3.0 - 0.5*y)))) &&
			       (h3 >=0.0 && h3<=0.5)) 
	      flag = 1;
	    /*region B*/   if ((((h1>=(1.0/3.0 - 0.5*y)&&h1<=(1.0/3.0-0.125*y)) &&
				 (h2>=0.5*h1 && h2<=(2.0*h1-(1.0/3.0) + 0.5*y))) ||
				((h1>=(1.0/3.0 - 0.125*y) && h1<=(1.0/3.0 + 0.25*y))
				 && (h2>=0.5*h1 && h2<=(1.0/3.0 + 0.25*y)))) && 
			       (h3 >= 0.0 && h3 <= (1.0/6.0 - (2.0/9.0)*y*(3*h1-1)))) 
	      flag = 1;
	    /*region D*/   if ((((h1>=(1.0/3.0-0.125*y) && h1<=1.0/3.0) && (h2>=(1.0/3.0
										 + 0.25*y) && h2 <= (2*h1-(1.0/3.0) + 0.5*y))) ||
				((h1>=1.0/3.0 && h1<=(1.0/3.0+0.25*y)) && (h2>=(1.0/3.0+
										0.25*y) && h2<=(2.0/3.0 + 0.5*y - h1)))) && (h3 >=
															     ((2.0/9.0)*y*(3*h2-1) - 1.0/6.0) && h3<=(1.0/6.0 -
																				      (2.0/9.0)*y*(3*h1-1))))   
	      flag = 1;
	    /*region E*/   if ((h2>=(1.0/3.0 - 0.5*y) && h2<=(1.0/3.0 + 0.25*y)) &&
			       (h1>=0.5*h2 && h1<=(0.5*h2 + 1.0/6.0 -0.25*y)) &&
			       (h3 >=0.0 && h3<= (1.0/3.0 - (1.0/9.0)*y*(3*h2-1))))  
	      flag = 1;
	    /*region F*/   if ((h2>=(1.0/3.0 + 0.25*y) && h2<=(1.0/3.0 + 0.5*y)) &&
			       (h1 >= 0.5*h2 && h1<=(0.5*h2 + 1.0/6.0 - 0.25*y)) &&
			       (h3>=((2.0/9.0)*y*(3*h2-1) - 1.0/6.0) && h3<=(1.0/3.0 -
									     (1.0/9.0)*y*(3*h2-1))))   
	      flag = 1;
	    if (flag == 1) {
	      w = 12;
	      if (h2==2*h1 || h2 == 0.5*h1) 
		w = 6;
	      if (h1 == h2 && (h3 == 0.0 || h3 == 0.5))  
		w = 6;
	      if (h1 == 1.0/6.0 && h2 == 1.0/3.0 && h3 == 1.0/3.0)  
		w = 3;
	      if (h1 == 1.0/3.0 && h2 == 1.0/6.0 && h3 == 1.0/6.0)  
		w = 3;
	      if (h1 ==0.0 && h2 == 0.0 )  
		w = 2;
	      if (h1 == 0.0 && h2 == 0.0 && (h3==0.0 || h3==0.5)) 
		w = 1;
		    fprintf(outfile,"%12.8f",r1);
		    fprintf(outfile,"%12.8f",r2);
		    fprintf(outfile,"%12.8f",r3);
		    fprintf(outfile,"%4.0d\n",w);
		 
	      num = num + 1;
	    }
	    r3 = r3 + del;
	  }
	  r2 = r2 + del;
	}
	r1 = r1 - del;
      }
    }      /*end of alpha < 90 section*/
    printf("Number of k points generated = %d\n",num);
  }         /*end of section where any calculation is done*/
  fclose(outfile);
}            /*end*/

