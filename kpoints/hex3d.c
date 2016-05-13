/* hex3d.c --  a program to generate special k point sets for 2-d hex */
/* systems.  Based on the Ramirez and Bohm method, 1986. */
/* written by John Mitchell, 2/2/89 */
/* Modified for EHMACC input file format,  6/23/89 -- J.M. 
   Modified to work with bind 03/19/97 -- Greg Landrum */



/* This program basically takes the framework from hex2d.c and fleshes */
/* it out by adding the third dimension.  It does this by computing the */
/* two dimensional grid as if the system were in the plane, and then */
/* adds on a k3 component which is determined by evenly dividing the */
/* c* axis into as many points as desired by the user.  In a matter of */
/* speaking, the program just takes the hex2d points and adds a component */
/* from the one dimensional k point generator, chain.c (c.f.) */
/* ------------------------------------------------------------------ */

#include <stdio.h>
#include <math.h>

#define MAX 34      /* maximum value of nd */
/* with this, max number of k points */
/* in each plane || k1-k2 plane is */
/* 94,125 */
/* Of course, since the 3rd dimension is */
/* you can have infinite number in this direction */

/* function to return true (1) if a positive integer is even */
int Even(x)
  int x;
  
{
  int i;
  while (x > 1)
    x = x - 2;
  if (x == 1)
    i = 0;
  else
    i = 1;
  return(i);
} /* end function Even */

/* ---------------------------------------------------------------- */

/* function to test whether a point is in the irreducible wedge or not */
/* Note:  this only tests in the k1-k2 plane; the third dimension is */
/*        automatically ok due to the way the main loop is set up to */
/*        keep k3 in the range 0 < k3 < 0.5 -- strict inequalities! */

int test(i,j)
  double i,j;
  
{
  int t;
  
  if (i>=0.50001)
    t = 0;
  else
    if (((i<=(0.333333333333)) &&
	 (j+i<=(2*i))) ||
	((i>(0.333333333333)) &&
	 ((i+(j/(2.0)))<=0.500001)))
      t = 1;
    else
      t = 0;
  return(t);
}
/* end test */

/* ----------------------------------------------------------------- */


main()
  
{
  int i,n;    /* subscript variables */
  int nd;     /* number of points along G to K, input by user, must be even */
  int num;    /* number of k points generated */
  int q;      /* number of k points parallel to the c* axis */
  double k[MAX][MAX][3];  /* k points.  i labels the "root" point along GK,
			     n labels which point along line parallel to GM
			     starting at i. Naturally, 3 components for last index */
  int w[MAX][MAX];            /* weight associated with k point */
  int check;          /* boolean - true if point in zone, false if out */
  double alpha,zeta;  /* intermediate numbers used in the calculation */
  double del;
  double K,k1,k2,k3;
  int progflag;
  FILE *outfile;
  
  /* open the output file for writing */
  outfile = fopen("kout","w");
  /* some helpful info */
  printf("\n\n\n");
  printf("Generator of 3d hexagonal k points.  Your direct lattice vectors \n");
  printf("must form an angle of 120 degrees or these points will be wrong!\n");
  printf("You will enter two seed numbers next. \n");
  printf("First, you will be asked for an EVEN number; this is used to\n");
  printf("Generate a two dimensional grid of k points.  You will then be\n");
  printf("asked for the number of k points you want parallel to the c* axis.\n");
  printf("This may be either odd or even.\n");
  printf("\n\n");
  printf("\n\n In general you will get q*(3/8)*(n^2 + 2*n) points, where\n");
  printf("n is your seed number, which must be even, and q is the\n");
  printf("number of points parallel to the c* axis\n");
  printf("\n\n");
  /* ask user for the number of points he wants along GK */
  /* and keep doing it until it is even */
  nd = 1;
  while (Even(nd) == 0)
      {
	printf("Enter an EVEN number as a seed \n");
	scanf("%d",&nd);
      }
  /* now find out how many parallel to the c* axis */
  printf("How many points parallel to the c* axis ? \n");
  scanf("%d",&q);
  /* compute the spacing along c* */
  del = 0.5/q;
  /* now for the guts of the program */
  num = 0;
  for (i=0;i<nd;++i)  /* loop over the nd "root" points */
      {
	/* alpha takes different values for odd and even i */
	if (Even(i))
	  alpha = (2.0 + 3.0*i)/(9.0*nd);
	else
	  alpha = (1.0 + 3.0*i)/(9.0*nd);
	n = 0;
	do
	    {
	      /* zeta takes different values for odd and even i */
	      if (Even(i))
		zeta = (4.0*n)/(3.0*nd - 3.0*i - 2.0);
	      else
		zeta = (4.0*n)/(3.0*nd - 3.0*i - 1.0);
	      k[i][n][0] = 0.5*zeta*(1-3.0*alpha) + alpha;
	      k[i][n][1] = alpha;
	      /* loop over the k3 component */
	      k3 = del/2.0;
	      while(k3<0.5)
		  {
		    /* check to see if the point lies in the boundary defined by the two */
		    /* dimensional Brillouin Zone.  Note that the boundary at the top */
		    /* of the zone is taken care of by the loop conditions on k3 -- that is */
		    /* k3 is always greater than zero and less than 1/2, so the k3 component */
		    /* will always lie in the Brillouin Zone */
		    check = test(k[i][n][0],k[i][n][1]);
		    if (check)   /* assign weights and write the point to file */
			{
			  k1 = k[i][n][0];
			  k2 = k[i][n][1];
			  w[i][n] = 24;                 /* general points */
			  if (fabs(k1-k2)<0.00000001)
			    w[i][n]=12;               /* on plane containing GK */;
			  if (k2<0.00000001)
			    w[i][n]=12;               /* on plane containing GM */
			  if ((fabs(k1 - (1.0/3.0)) >= 0.000000001) &&
			      (fabs(k1+(k2/2.0)-0.5)<0.0000001))
			    w[i][n]=12;               /* on plane containing MK */
			  
			  /* Note that this k point set will never give any high symmetry points */
			  /* so I have not considered them in this weighting scheme */
			  /* Additionally, since 0 < k3 < 1/2 (strict inequalities), no points will */
			  /* lie on these mirror planes */
			  
			  /* now print the k points out to the file 'kout' */
			  fprintf(outfile,"%12.8f",k1);
			  fprintf(outfile,"%12.8f",k2);
			  fprintf(outfile,"%12.8f",k3);
			  fprintf(outfile,"%4.0d\n",w[i][n]);
			  
			  ++num;   /* since in the zone, increment total # of points by 1 */
			}
		    k3 += del;
		  }  /* end loop over k3 */
	      ++n;
	    }
	while (check != 0);  /* end loop over n */
      } /* end loop over i */
  printf ("Number of k points generated = %d \n",num);
  fclose(outfile);
}  /* end of program */


