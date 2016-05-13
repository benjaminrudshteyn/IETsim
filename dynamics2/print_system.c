/*******************************************************
 *      Copyright (C) 1995 Greg Landrum
 *
 *  This file is part of yaehmop.
 *
 *   This is free software.
 * 
 *  Permission is granted to modify, or otherwise fold, spindle, and mutilate this
 *    code provided all copyright notices are left intact.
 *
 *  This code may be distributed to your heart's content, in whatever form,
*    provided no fee is charged for the distribution, all copyright notices are
*    left intact, and the source is distributed (without fee) along with any
*    binaries to anyone who requests it.
*
*  There are, of course, no warranties at all on this program.
*
********************************************************************/



/****************************************************************************
*
*  created:  Robert Snoeberger  December 2008
*
*****************************************************************************/
#include "bind.h"


/****
  Recent edit history

****/


/****************************************************************************
 *
 *                   Procedure print_system 
 *
 * Arguments:  cell: pointer to cell type
 *          details: pointer to detail type
 *          overlapR: hermetian_matrix_type
 *            hamilR: hermetian_matrix_type
 *          overlapK: hermetian_matrix_type
 *            hamilK: hermetian_matrix_type
 *   cmplx_hamil, cmplx_overlap: pointers to complex
 *          eigenset: eigenset_type
 *        properties: pointer to prop_type
 *     avg_prop_info: pointer to avg_prop_info_type
 * work1,work2,work3: pointers to reals
 *        cmplx_work: pointer to complex
 *          num_orbs: int
 * orbital_lookup_table: pointer to int.
 *
 * Returns: none
 *
 *
 *
 *
 ****************************************************************************/
void print_system(cell,details,overlapR,hamilR,overlapK,hamilK,
			cmplx_hamil,cmplx_overlap,
			eigenset,work1,work2,work3,cmplx_work,
			properties,avg_prop_info,
			num_orbs,orbital_lookup_table,file_name)
  cell_type *cell;
  detail_type *details;
  hermetian_matrix_type overlapR,hamilR;
  hermetian_matrix_type overlapK,hamilK;
  complex *cmplx_hamil,*cmplx_overlap;
  eigenset_type eigenset;
  real *work1,*work2,*work3;
  complex *cmplx_work;
  prop_type *properties;
  avg_prop_info_type *avg_prop_info;
  int num_orbs;
  int *orbital_lookup_table;
  char *file_name;
{
  int i,j,k,l,m;
  FILE *f;
  

  /* Print Hamiltonian to file hamiltonian.dat */
  f = fopen("hamiltonian.dat","w");
  fprintf(f,"%d\n",num_orbs);
  k = 0;
  for (i=0;i<num_orbs;i++){
    for (j=0;j<num_orbs;j++){
      fprintf(f,"%22.15e ",HERMETIAN_R(hamilR,i,j));
      k += 1;
      if (k==5) {
	 fprintf(f,"\n");
         k = 0;
      }	 
    }
    fprintf(f,"\n");
    k = 0;
  }
  fclose(f);


  /* Print Overlap to file overlap.dat */
  f = fopen("overlap.dat","w");
  fprintf(f,"%d\n",num_orbs);
  k = 0;
  for (i=0;i<num_orbs;i++){
    for (j=0;j<num_orbs;j++){
      fprintf(f,"%22.15e ",HERMETIAN_R(overlapR,i,j));
      k += 1;
      if (k==5) {
	 fprintf(f,"\n");
         k = 0;
      }	 
    }
    fprintf(f,"\n");
    k = 0;
  }
  fclose(f);

  if (!details->just_matrices) { /* If System was diagonalized */
    /* Print Orbital to file orbital.dat */
    f = fopen("orbital.dat","w");
    fprintf(f,"%d\n",num_orbs);
    k = 0;
    for (i=0;i<num_orbs;i++){
      for (j=0;j<num_orbs;j++){
        fprintf(f,"%22.15e ",EIGENVECT_R(eigenset,j,i));
        k += 1;
        if (k==5) {
      	   fprintf(f,"\n");
           k = 0;
        }	 
      }
      fprintf(f,"\n");
      k = 0;
    }
    fclose(f);


    
    /* Print energies to file energy.dat */
    f = fopen("energy.dat","w");
    fprintf(f,"%d\n",num_orbs);
    for (i=0;i<num_orbs;i++){
      fprintf(f,"%22.15e\n",eigenset.val[i]);
    }
    fclose(f);
    
  }


} 



