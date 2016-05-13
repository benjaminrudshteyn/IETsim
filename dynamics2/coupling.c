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
*  created:  Robert Snoeberger  January 2008
*
*****************************************************************************/
#include "bind.h"

real gauss(real E,real Ei,real width);

/****
  Recent edit history

****/


/****************************************************************************
 *
 *                   Procedure dynamics
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
void ecoupling(cell,details,overlapR,hamilR,overlapK,hamilK,
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
  static char tempfilename[240];
  static char *label=0;
  static FILE *sparse_OVfile,*sparse_HAMfile;
  k_point_type *kpoint;
  real *mat_save;
  int i,j,k,l,m;
  int itab,jtab,ktab;
  int ltab,mtab;
  int diag_error;
  real temp;
  real total_energy,tot_chg;
  int electrons_so_far;
  real *occupations;
  real *chg_mat;
  int num_KPOINTS;
  int overlap_file,hamil_file;
#ifdef USE_LAPACK
  char jobz, uplo;
  int info, itype;
  int num_orbs2;
#endif
  hermetian_matrix_type adsorb_over,adsorb_ham; /* overlap and hamiltonian matrix for adsorbate */
  hermetian_matrix_type acceptor_over,acceptor_ham; /* overlap and hamiltonian matrix for acceptor */
  hermetian_matrix_type big_over,big_ham,big_S; /* overlap, hamiltonian and vector matrix for diabatic state */
  hermetian_matrix_type one,two,three; 
  eigenset_type adeigenset,acceptorset;
  BOOLEAN ad1,ad2;
  int num_adorbs,num_acceptorbs;
  real twopiohbar = 9.5458397443104435; // units: 1/(eV * fs)
  real rate,width,Ei;
  FILE *f;
  
  fprintf(output_file,"\n# ************Electronic Coupling*****************\n");
 
  /* IMPORTANT!!!!!!!!!!!!!!!!!!!!
 *   
 *   To make programming easier, it is assummed that the orbitals are in
 *   the order of adsorbate orbitals first and then acceptor obitals second
 *   so the adsorbate atoms need to be the first atoms specified in the
 *   input file.  */

 
  /* 
   * Determine number of orbitals in adsorbate and acceptor system
   * */

  fprintf(output_file,"\nThere are %d adsorbate atoms and %d acceptor atoms\n",
                      details->adsorbate_num_atoms,
                      cell->num_atoms - details->adsorbate_num_atoms);
  
  fprintf(output_file,"\nAdsorbate Atom Table\n");
  /* determine number of adsorbate orbitals*/
  details->adsorbate_num_orbs = 0;
  for (i=0;i<details->adsorbate_num_atoms;i++){
     fprintf(output_file,"%3d,  ",details->adsorbate[i]);
     if (i%13==12) fprintf(output_file,"\n");
     k = details->adsorbate[i];
     //details->adsorbate_num_orbs += cell->atoms[i].num_valence;
     if (k+1>cell->num_atoms-1) {
       details->adsorbate_num_orbs += hamilR.dim - orbital_lookup_table[k];
     } else {
       details->adsorbate_num_orbs += orbital_lookup_table[k+1]-orbital_lookup_table[k];
     }
  }
  
  num_adorbs = details->adsorbate_num_orbs;
  num_acceptorbs = num_orbs - details->adsorbate_num_orbs;

  fprintf(output_file,"\n\nThere are %d adsorbate orbitals and %d acceptor orbitals\n",num_adorbs,num_acceptorbs);




  /* allocate memory for adsorbate orbital table */
  details->adsorbate_orbs = (int *)calloc(details->adsorbate_num_orbs,sizeof(int));

  /* Allocate memory for Adsorbate and Acceptor Matricies */
  adsorb_over.dim = details->adsorbate_num_orbs;
  adsorb_over.mat = (real *)calloc(adsorb_over.dim*adsorb_over.dim,sizeof(real));

  adsorb_ham.dim = details->adsorbate_num_orbs;
  adsorb_ham.mat = (real *)calloc(adsorb_ham.dim*adsorb_ham.dim,sizeof(real));

  adeigenset.vectR = (real *)calloc(num_adorbs*num_adorbs,sizeof(real));
  adeigenset.val = (real *)calloc(num_adorbs,sizeof(real));

  acceptor_over.dim = num_acceptorbs;
  acceptor_over.mat = (real *)calloc(acceptor_over.dim*acceptor_over.dim,
                        sizeof(real));
  
  acceptor_ham.dim = num_acceptorbs;
  acceptor_ham.mat = (real *)calloc(acceptor_ham.dim*acceptor_ham.dim,sizeof(real));

  acceptorset.vectR = (real *)calloc(num_acceptorbs*num_acceptorbs,sizeof(real));
  acceptorset.val = (real *)calloc(num_acceptorbs*num_acceptorbs,sizeof(real));  

  big_over.dim = num_orbs;
  big_over.mat = (real *)calloc(big_over.dim*big_over.dim,sizeof(real));
 
  big_ham.dim = num_orbs;
  big_ham.mat = (real *)calloc(big_ham.dim*big_ham.dim,sizeof(real));

  big_S.dim = num_orbs;
  big_S.mat = (real *)calloc(big_S.dim*big_S.dim,sizeof(real));




  /* Make the adsorbate orbital table */
  fprintf(output_file,"\nAdsorbate Orbital Table\n");
  m = 0;
  for (i=0;i<details->adsorbate_num_atoms;i++){
     k = details->adsorbate[i];
     if (k+1>cell->num_atoms-1) {
       l = num_orbs - orbital_lookup_table[k];
     } else {
       l = orbital_lookup_table[k+1]-orbital_lookup_table[k];
     }
     for (j=0;j<l;j++){
       details->adsorbate_orbs[m] = orbital_lookup_table[k]+j;
       fprintf(output_file,"%3d,  ",details->adsorbate_orbs[m]);
       m += 1;
       if (m%13==12){
	 fprintf(output_file,"\n");
       }
     }
  }


  
  /* 
   * Make adsorbate Hamiltonian and Overlap 
   * */	 

  /* Old way which ignores orbital order */
  //for (i=0;i<num_adorbs;i++){
  //  for (j=i;j<num_adorbs;j++){
  //    k = details->adsorbate_orbs[i];
  //    l = details->adsorbate_orbs[j];
  //    //printf("Ha[%d,%d] = H[%d,%d] = %f\n",i,j,k,l,MATRIX(hamilR,k,l));
  //    MATRIX(adsorb_ham,i,j) = MATRIX(hamilR,k,l);
  //    MATRIX(adsorb_over,i,j) = MATRIX(overlapR,k,l);
  //  }
  //}

  /* Easy way which assumes adsorbate orbitals are the first #(num_adorbs) orbitals */
  for (i=0;i<num_adorbs;i++){
    for (j=i;j<num_adorbs;j++){
      MATRIX(adsorb_ham,i,j) = MATRIX(hamilR,i,j);
      MATRIX(adsorb_over,i,j) = MATRIX(overlapR,i,j);

      MATRIX(big_over,i,j) = MATRIX(overlapR,i,j); /* Overlap Matrix for Diabat */
      MATRIX(big_over,j,i) = MATRIX(overlapR,i,j);
      /* Needed since we need to perform S*Over*H*S^T later */ 
    }
  }


  for (i=0;i<num_acceptorbs;i++){
    for (j=i;j<num_acceptorbs;j++){
      k = num_adorbs+i;
      l = num_adorbs+j;
      MATRIX(acceptor_ham,i,j) = MATRIX(hamilR,k,l);
      MATRIX(acceptor_over,i,j) = MATRIX(overlapR,k,l);

      MATRIX(big_over,k,l) = MATRIX(overlapR,k,l);
      MATRIX(big_over,l,k) = MATRIX(overlapR,k,l);
    }
  }

  /*
   * Diagonalize Adsorbate and Acceptor System
   */
      /*******
	
	now diagonalize that beast by calling the FORTRAN subroutine used
	to diagonalize stuff in new3 and CACAO.
	
	THIS REALLY SHOULD BE REPLACED with a routine written in C, so if you
	happen to have some time on your hands....
	
	********/
  printf("}");
  cboris(&(num_adorbs),&(num_adorbs),adsorb_ham.mat,adsorb_over.mat,adeigenset.vectR,adeigenset.val,work1,
	     work2,&diag_error);
  printf("{\n");

  fprintf(output_file,"\nAdsorbate Energies  %d\n",diag_error);
  for (i=0;i<num_adorbs;i++){
    fprintf(output_file,"%d:---> %f\n",i+1,adeigenset.val[i]);
  }
  

  printf("}");
  cboris(&(num_acceptorbs),&(num_acceptorbs),acceptor_ham.mat,acceptor_over.mat,acceptorset.vectR,acceptorset.val,
              work1,work2,&diag_error);
  printf("{\n");

  fprintf(output_file,"\nAcceptor Energies  %d\n",diag_error);
  for (i=0;i<num_acceptorbs;i++){
    fprintf(output_file,"%d:---> %f\n",i+1,acceptorset.val[i]);
  }

  /* Diagonalization should be complete - If so the hamiltoniam matricies will contain the eigenvectors*/
  /*
  printf("Checking Adsorbate Eigenvectors \n");
  for (i=0;i<num_adorbs;i++){
      printf("%f  %f\n",MATRIX(adsorb_ham,i,0),MATRIX(adsorb_ham,i,1));
  }

  printf("Checking Acceptor Eigenvectors \n");
  for (i=0;i<num_acceptorbs;i++){
      printf("%f  %f\n",MATRIX(acceptor_ham,i,0),MATRIX(acceptor_ham,i,1));
  }
  */


  /* Need to make matrix of diabat eigenvectors */
  
  for (i=0;i<num_adorbs;i++){
    for (j=0;j<num_adorbs;j++){
       MATRIX(big_S,i,j)  = MATRIX(adsorb_ham,i,j);
    } 
  }
 
  for (i=0;i<num_acceptorbs;i++){
    for (j=0;j<num_acceptorbs;j++){
      k = num_adorbs+i;
      l = num_adorbs+j;
      MATRIX(big_S,k,l) = MATRIX(acceptor_ham,i,j);
    }
  }
  
  /* copy hamiltonian matrix to big_ham so it can be saved for later if needed */
  for (i=0;i<num_orbs;i++){
    for (j=i;j<num_orbs;j++){
      MATRIX(big_ham,i,j) = MATRIX(hamilR,i,j);
      MATRIX(big_ham,j,i) = MATRIX(hamilR,i,j);
    }
  }
  /*
  printf("Checking Hamiltonian \n");
  for (i=0;i<num_orbs;i++){
      printf("%f  %f %f  %f\n",MATRIX(big_ham,i,0),MATRIX(big_ham,i,1),MATRIX(big_ham,i,2),MATRIX(big_ham,i,3));
  }

  printf("Checking Diabat Overlap \n");
  for (i=0;i<num_orbs;i++){
      printf("%f  %f %f  %f\n",MATRIX(big_over,i,0),MATRIX(big_over,i,1),MATRIX(big_over,i,2),MATRIX(big_over,i,3));
  }

  printf("Checking Diabat Eigenvectors \n");
  for (i=0;i<num_orbs;i++){
      printf("%f  %f %f  %f\n",MATRIX(big_S,i,0),MATRIX(big_S,i,1),MATRIX(big_S,i,2),MATRIX(big_S,i,3));
  }
 */

  /* Need to perform Algebra T*H*T^  */

  /* Need working space */
  one.dim = num_orbs;
  one.mat = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  two.dim = num_orbs;
  two.mat = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  three.dim = num_orbs;
  three.mat = (real *)calloc(num_orbs*num_orbs,sizeof(real));


  /* one = big_ham * big_S^ */
  for (i=0;i<num_orbs;i++){
    for (j=0;j<num_orbs;j++){
      MATRIX(one,i,j) = 0.0;
      for (k=0;k<num_orbs;k++){
        MATRIX(one,i,j) += MATRIX(big_ham,i,k)*MATRIX(big_S,j,k);
      }
    }
  }

  /* three = big_S * two */
  for (i=0;i<num_orbs;i++){
    for (j=0;j<num_orbs;j++){
      MATRIX(three,i,j) = 0.0;
      for (k=0;k<num_orbs;k++){
        MATRIX(three,i,j) += MATRIX(big_S,i,k)*MATRIX(one,k,j);
      }
    }
  }  

  /* three should be diagonalized to diabatic state */
  
  fprintf(output_file,"\nDiabat Energies\n");
  for (i=0;i<num_orbs;i++){
    fprintf(output_file,"%d:---> %f\n",i+1,MATRIX(three,i,i));
  }


  /* Determine the coupling to the Adsorbate State */
  fprintf(output_file,"\nElectronic Coupling to Adsorbate State %d\n",details->coupling_state);
  i = details->coupling_state-1;
  for (j=num_adorbs;j<num_orbs;j++){
     fprintf(output_file,"%d C--> %-9.6f   %-9.6f\n",j-num_adorbs+1,pow(MATRIX(three,i,j),2.0),MATRIX(three,j,j));
  }

  /* Calculate Rate using Fermi's Golden Rule */
  width = 0.693147181/pow((details->coupling_width)/2.0,2.0); 
  i = details->coupling_state-1;
  Ei = MATRIX(three,i,i);
  rate = 0.0;
  for (j=num_adorbs;j<num_orbs;j++){
    rate += twopiohbar*pow(MATRIX(three,i,j),2.0)*gauss(MATRIX(three,j,j),
                          Ei,width);
  } 
  printf("\n\nThe IET rate is %f fs\n",1/rate);
} 



/******************* Gaussian Function ****************/
real gauss(real E,real Ei,real width){
  return exp(-width*pow(E-Ei,2.0));
}

