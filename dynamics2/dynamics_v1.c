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
void propagate(eigenset_type a,wavepacket_type b,real c);
real survival(wavepacket_type a,wavepacket_type b);
void wave_header(char *filename,cell_type *cell,int i,int *j);
void wave_packet(char *filename,wavepacket_type ao);

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
 * Action:  Electron Dynamics 
 *     things are done:
 *    - write any required info to the output file
 *
 *   Some helpful definitions:
 *
 *    
 *              N                     |              N                     
 *            -----                   |            -----                   
 *             \                      |             \                      
 *    H(k) :=   )   exp(-i k . r) H(r)|    S(k) :=   )   exp(-i k . r) S(r)   
 *             /                      |             /                      
 *            -----                   |            -----                   
 *            r = 1                   |            r = 1                   
 *                                    |
 *                                    |
 *
 *   The work arrays are used as temporary memory in the various functions called
 *    by this one.  The dimensions should be:
 *         work1,work2: num_orbs;
 *    work3,cmplx_work: num_orbs*num_orbs;
 *
 *
 *
 ****************************************************************************/
void dynamics(cell,details,overlapR,hamilR,overlapK,hamilK,
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
  wavepacket_type ao,mo,save_mo;  /* atomic orbital and molecular orbital for wavepacket propagation */
  hermetian_matrix_type adsorb_over,adsorb_ham; /* overlap and hamiltonian matrix for adsorbate */
  eigenset_type adeigenset;
  real norm,normad;  /* For normalizing eigenvectors */
  real time;         /* Time during dynamics */
  BOOLEAN ad1,ad2;
  int num_adorbs;
  int istep;
  static char cubefile[20],wavefile[20];
  FILE *f;
  
  fprintf(output_file,"\n# ************ELECTRON DYNAMICS*******************\n");
  
  /* 
   * Allocate memory for wavepacket  
   * */  
  ao.dim = hamilR.dim;
  ao.val = (complex *)calloc(ao.dim,sizeof(complex));
  mo.dim = hamilR.dim;
  mo.val = (complex *)calloc(mo.dim,sizeof(complex));
  save_mo.dim = mo.dim;
  save_mo.val = (complex *)calloc(mo.dim,sizeof(complex));


  /*-----------------Prepare Initial Wavepacket---------------------------*/
  
  /* 
   * Create adsorbate orbital table
   * */

  fprintf(output_file,"\nThere are %d adsorbate atoms\n",details->adsorbate_num_atoms);
  
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
  fprintf(output_file,"\n\nThere are %d adsorbate orbitals\n",details->adsorbate_num_orbs);
  
  /* allocate memory for adsorbate orbital table and adsorbate system */
  num_adorbs = details->adsorbate_num_orbs;
  details->adsorbate_orbs = (int *)calloc(details->adsorbate_num_orbs,sizeof(int));
  adsorb_over.dim = details->adsorbate_num_orbs;
  adsorb_over.mat = (real *)calloc(adsorb_over.dim*adsorb_over.dim,sizeof(real));
  adsorb_ham.dim = details->adsorbate_num_orbs;
  adsorb_ham.mat = (real *)calloc(adsorb_ham.dim*adsorb_ham.dim,sizeof(real));
  adeigenset.vectR = (real *)calloc(num_adorbs*num_adorbs,sizeof(real));
  adeigenset.val = (real *)calloc(num_adorbs,sizeof(real));

  /* finally make the adsorbate orbital table */
  fprintf(output_file,"\nAdsorbate Orbital Table\n");
  m = 0;
  for (i=0;i<details->adsorbate_num_atoms;i++){
     k = details->adsorbate[i];
     if (k+1>cell->num_atoms-1) {
       l = hamilR.dim - orbital_lookup_table[k];
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
  for (i=0;i<num_adorbs;i++){
    for (j=i;j<num_adorbs;j++){
      k = details->adsorbate_orbs[i];
      l = details->adsorbate_orbs[j];
      //printf("Ha[%d,%d] = H[%d,%d] = %f\n",i,j,k,l,MATRIX(hamilR,k,l));
      MATRIX(adsorb_ham,i,j) = MATRIX(hamilR,k,l);
      MATRIX(adsorb_over,i,j) = MATRIX(overlapR,k,l);
    }
  }
  
  //printf("\nChecking Overlap\n");
  //for (i=0;i<num_adorbs;i++){
  //  for (j=0;j<num_adorbs;j++){
  //    printf("%f   ",MATRIX(adsorb_over,i,j));
  //  }
  //  printf("\n");
  //}
    
  //printf("\nChecking Hamil\n");
  //for (i=0;i<num_adorbs;i++){
  //  for (j=0;j<num_adorbs;j++){
  //    printf("%f   ",MATRIX(adsorb_ham,i,j));
  //  }
  //  printf("\n");
  //}

  /*
   * Diagonalize Adsorbate System
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

  /* Move eigenvectors to adeigenset */
  //for (i=0;i<hamilR.dim;i++){
  //  for (j=0;j<hamilR.dim;j++){
  //    TRANSFORM(adeigenset,i,j) = MATRIX(adsorb_ham,i,j);
  //  }
  //}
  
  fprintf(output_file,"\nAdsorbate Energies\n");
  for (i=0;i<num_adorbs;i++){
    fprintf(output_file,"%d:---> %f\n",i+1,adeigenset.val[i]);
  }
  
  /* 
   *  Normalize eigenvectors (complete system and adsorbate) to T**H*T = I
   * */
  for (i=0;i<num_orbs;i++){
    norm = 0.0;
    for (j=0;j<num_orbs;j++){
      norm +=  TRANSFORM(eigenset,i,j)*TRANSFORM(eigenset,i,j);
    }
    norm = pow(norm,-0.5);
    for (j=0;j<hamilR.dim;j++){
      TRANSFORM(eigenset,i,j) *= norm;
    }
  }

  for (i=0;i<num_adorbs;i++){
    normad = 0.0;
    for (j=0;j<num_adorbs;j++){
      normad += MATRIX(adsorb_ham,i,j)*MATRIX(adsorb_ham,i,j);
    }
    normad = pow(normad,-0.5);
    for (j=0;j<num_adorbs;j++){
      MATRIX(adsorb_ham,i,j) *= normad;
    }
  }
  
  /* Check eigenvectors */
  //printf("\nSystem eigenvectors\n");
  //for (i=0;i<num_orbs;i++){
  //  for (j=0;j<num_orbs;j++){
  //    printf("%f   ",TRANSFORM(eigenset,i,j));
  //  }
  //  printf("\n");
  //}
  //printf("\nAdsorbate eigenvectors\n");
  //for (i=0;i<num_adorbs;i++){
  //  for (j=0;j<num_adorbs;j++){
  //    //printf("%f   ",TRANSFORM(adeigenset,i,j));
  //    printf("%f   ",MATRIX(adsorb_ham,i,j));
  //  }
  //  printf("\n");
  //}

  /*
   * Initalize Wavepacket *
   * */
  k = details->initial_orbital-1;
  for (i=0;i<ao.dim;i++){
    ao.val[i]=0.0;
    mo.val[i]=0.0;
    save_mo.val[i]=0.0;
  }
  for (i=0;i<num_adorbs;i++){
    //printf("adsorb orb %d\n",details->adsorbate_orbs[i]);
    ao.val[details->adsorbate_orbs[i]] = MATRIX(adsorb_ham,k,i);
  }
  //printf("\nCheck initial wavepacket\n");
  //for (i=0;i<ao.dim;i++){
  //  printf("%f\n",creal(ao.val[i]));
  //} 

  /* Transform Initial wavepacet from ao basis -> mo basis*/
  norm = 0.0;
  for (i=0;i<num_orbs;i++){
    for (j=0;j<num_orbs;j++){
      mo.val[i] += TRANSFORM(eigenset,i,j)*ao.val[j];
    }
    norm += creal(mo.val[i]*mo.val[i]);
  }        
  norm = pow(norm,-0.5);   /* Make sure mo is normalized to one */
  for (i=0;i<num_orbs;i++){ /* It is strange that the mo isnt normalized - is diagonalization to inaccurate? */
    mo.val[i] *= norm;
    save_mo.val[i] = mo.val[i];
  }

  //printf("\nCheck initial mo wavepacket\n");
  //for (i=0;i<mo.dim;i++){
  //  printf("%f   %f\n",creal(mo.val[i]),cimag(mo.val[i]));
  //} 

  /* 
   * The system is now ready for dynamics!!!
   * */

  /*--------------Electron Dynamics----------------------------------*/
  f = fopen(file_name,"w");
  fprintf(f,"# Electron Dynamics\n# Time  Survival\n");
  
  istep = 0;
  time = 0.0;
  /* Write Initial Wavepacket Data to file */
  fprintf(f,"%f   %f\n",time,survival(mo,save_mo));

  /* prepare wave filename */
  sprintf(wavefile,"%s.wave",file_name);
  printf("The wave file is named %s\n",&wavefile);

  /* call function to write header for wavefile */
  wave_header(wavefile,cell,num_orbs,orbital_lookup_table);
  /* write wavepacket to wave file */
  wave_packet(wavefile,ao);
  
  if (details->cube){
    /* prepare cube filename */
    sprintf(cubefile,"%s.0.cube",file_name,istep);
    printf("The cube file is named %s\n",&cubefile);
  
    /* call function to generate cube */
    write_cube(cubefile,ao,overlapR,cell,orbital_lookup_table,details);
  }
    
  while (time < details->time_length){
     
    /*
     *  Propagate mo 
     *  */
    printf("|t=%3.2f> = ",time+details->timestep);
    propagate(eigenset,mo,details->timestep);
    printf("exp^(-iHt/hbar)*|t=%3.2f>\n",time); 
    time += details->timestep;
    istep += 1;
    
    /* Write Time-evolved Wavepacket Data to file */
    fprintf(f,"%f   %f\n",time,survival(mo,save_mo));
    
	   
    /* Transform wavepacket from mo basis -> ao basis */
    for (i=0;i<num_orbs;i++){
      ao.val[i] = 0.0;
      for (j=0;j<num_orbs;j++){
        ao.val[i] += TRANSFORM(eigenset,j,i)*mo.val[j];
      }
    }
    /* write wavepacket to wave file */
    wave_packet(wavefile,ao);
    
    /* 
     * Write out cube if desired
     * */
    if (details->cube){

      /* prepare cube filename */
      sprintf(cubefile,"%s.%d.cube",file_name,istep);
      printf("The cube file is named %s\n",&cubefile);
      
      /* call function to generate cube */
      write_cube(cubefile,ao,overlapR,cell,orbital_lookup_table,details);
    }
    /* Time loop finished */
  } 
  
  /*-------------Done with Dynamics----------------------------------*/
  fclose(f);
}
  

void propagate(eigen,psi,timestep)
  eigenset_type eigen;
  wavepacket_type psi;
  real timestep;
{
  real hbar=0.65821189916; // eV*fs (electronVolt femtosecond)
  real save_psi;
  int i,j,k;

  /* 
   * Apply Time evolution operator
   * cos(E*tau/hbar)-isin(E*tau/hbar)
   * */
  for (i=0;i<psi.dim;i++){
    psi.val[i] = psi.val[i]*cexp(-I*ENERGY(eigen,i)*timestep*(1/hbar));
  }
  /* The wavepacket has evolved! */

}

real survival(mo,save_mo)
  wavepacket_type mo,save_mo;
{
  complex surv;
  int i,j;

  /*
   * Calculate <t=0|t=time> 
   * */
  surv=0.0;
  for (i=0;i<mo.dim;i++){
    surv += conj(save_mo.val[i])*mo.val[i];
  }
  return cabs(surv);
}

void wave_header(filename,cell,num_orbs,orbital_lookup_table)
  char *filename;
  cell_type *cell;
  int num_orbs;
  int *orbital_lookup_table;
{
  int i,j,k;
  int iorbs;  
  FILE *f;

  /* Open wave file */
  f = fopen(filename,"w");
  
  /* Write Molecule information to wavefile*/
  fprintf(f,"MOLECULE\n%d\n",cell->num_atoms);
  for (i=0;i<cell->num_atoms;i++){
    fprintf(f,"%3d  %12.8f  %12.8f  %12.8f\n",cell->atoms[i].at_number,cell->atoms[i].loc.x,
		    cell->atoms[i].loc.y,cell->atoms[i].loc.z); 
  }
  fprintf(f,"\n");

  /* Write BASIS information to wavefile */
  fprintf(f,"BASIS\n%d\n",num_orbs);
  for (i=0;i<cell->num_atoms;i++){
    if (i+1==cell->num_atoms) {
      iorbs = num_orbs - orbital_lookup_table[i];
    } else {
      iorbs = orbital_lookup_table[i+1] - orbital_lookup_table[i];
    }
    for (j=1;j<=iorbs;j++){
      /* atom n ml coef1 exp1 coef2 exp2  */
      if (j==1){	    
        fprintf(f,"%d  %d  %d  %f  %f  %f  %f\n",i,cell->atoms[i].ns,j,1.0,cell->atoms[i].exp_s,0.0,0.0);
      } else if (j<5){
        fprintf(f,"%d  %d  %d  %f  %f  %f  %f\n",i,cell->atoms[i].np,j,1.0,cell->atoms[i].exp_p,0.0,0.0);
      } else if (j<10){
        fprintf(f,"%d  %d  %d  %f  %f  %f  %f\n",i,cell->atoms[i].nd,j,cell->atoms[i].coeff_d1,
			cell->atoms[i].exp_d,cell->atoms[i].coeff_d2,cell->atoms[i].exp_d2);
      }	else {      
        fprintf(f,"%d  %d  %d  %f  %f  %f  %f\n",i,cell->atoms[i].nf,j,cell->atoms[i].coeff_f1,
			cell->atoms[i].exp_f,cell->atoms[i].coeff_f2,cell->atoms[i].exp_f2);
      }
    }
  }
  fprintf(f,"\n");
  
  /* Close wavefile */
  fclose(f);
}  

void wave_packet(filename,ao)
  char *filename;
  wavepacket_type ao;
{
  int i;
  FILE *f;

  /* Open wave file */
  f = fopen(filename,"a");

  /* Write WAVEPACKET information to wavefile */
  fprintf(f,"WAVEPACKET\n%d\n",ao.dim);
  for (i=0;i<ao.dim;i++){
    fprintf(f,"%16.12lf  %16.12lf\n",creal(ao.val[i]),cimag(ao.val[i]));
  }
  fprintf(f,"\n");

  /* Close wavefile */
  fclose(f);
}


