/* hex2d.c --  a program to generate special k point sets for 2-d hex */
/* systems.  Based on the Ramirez and Bohm method, 1986. */
/* written by John Mitchell, 2/1/89 */
/* Modified for EHMACC input file format, 6/23/89 -- J.M. 
   Modified to work with bind 03/19/97 -- Greg Landrum 
   */
/* ------------------------------------------------------------------ */

#include <stdio.h>
#include <math.h>

#define MAX 40      /* maximum value of nd */
/* with this, max number of k points = 94,125 */

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
  double k[MAX][MAX][2];  /* k points.  i labels the "root" point along GK,
			     n labels which point along line parallel to GM
			     starting at i. Naturally, 2 components for last index */
  int w[MAX][MAX];            /* weight associated with k point */
  int check;          /* boolean - true if point in zone, false if out */
  double alpha,zeta;  /* intermediate numbers used in the calculation */
  double k1,k2,k3;
  int progflag;
  FILE *outfile;
  
  /* open the output file for writing */
  outfile = fopen("kout","w"); 
  /* some helpful info */
  printf("\n\n\n");
  printf("Generator of 2d hexagonal k points.  Your direct lattice vectors \n");
  printf("must form an angle of 120 degrees or these points will be wrong!\n");
  printf("You will enter a seed number next.  Here are some helpful hints: \n");
  printf("\n\n");
  printf("         SEED            NUMBER OF K POINTS\n\n");
  printf("          2                       3\n");
  printf("          4                       9\n");
  printf("          6                      18\n");
  printf("          8                      30\n");
  printf("         10                      45\n");
  printf("         12                      63\n");
  printf("         14                      84\n");
  printf("\n\n In general you will get (3/8)*(n^2 + 2*n) points, where\n");
  printf("n is your seed number, which must be even.\n");
  printf("\n\n");
  /* ask user for the number of points he wants along GK */
  /* and keep doing it until it is even */
  nd = 1;
  while (Even(nd) == 0)
      {
	printf("Enter an EVEN number as a seed \n");
	scanf("%d",&nd);
      }
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
	      check = test(k[i][n][0],k[i][n][1]);
	      if (check)   /* assign weights and write the point to file */
		  {
		    k1 = k[i][n][0];
		    k2 = k[i][n][1];
		    k3 = 0.0;
		    w[i][n] = 12;                 /* general points */
		    if (fabs(k1-k2)<0.00000001)
		      w[i][n]=6;               /* on GK */;
		    if (k2<0.00000001)
		      w[i][n]=6;               /* on GM */
		    if ((fabs(k1 - (1.0/3.0)) >= 0.000000001) &&
			(fabs(k1+(k2/2.0)-0.5)<0.0000001))
		      w[i][n]=6;               /* on MK */
		    
		    /* Note that this k point set will never give G, M, or K -- so I have not
		       considered them in this weighting scheme */
		    
		    /* now print the k points out to the file 'kout' */
		    fprintf(outfile,"%12.8f",k1);
		    fprintf(outfile,"%12.8f",k2);
		    fprintf(outfile,"%12.8f",k3);
		    fprintf(outfile,"%4.0d\n",w[i][n]);
		    ++num;   /* since in the zone, increment total # of points by 1 */
		  }
	      ++n;
	    }
	while (check != 0);
      } /* end loop over i */
  printf ("Number of k points generated = %d \n",num);
  fclose(outfile);
}  /* end of program */
