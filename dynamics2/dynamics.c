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
void propagate(eigenset_type eigen,detail_type *details,wavepacket_type *ao,wavepacket_type mo,real time);
real survival(wavepacket_type a,hermetian_matrix_type overlapR,detail_type *details);
real normalization(wavepacket_type ao,hermetian_matrix_type overlapR);
real occupation(wavepacket_type ao,hermetian_matrix_type overlapR,detail_type *details);
void wave_header(char *filename,cell_type *cell,int i,int *j);
void wave_packet(char *filename,wavepacket_type ao);
void write_restart(char *restartfile,wavepacket_type ao,wavepacket_type mo,real *energy);
void project(wavepacket_type *projection,wavepacket_type ao,eigenset_type eigenset,hermetian_matrix_type overlapR);

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
  wavepacket_type ao,save_ao,mo,save_mo;  /* atomic orbital and molecular orbital for wavepacket propagation */
  hermetian_matrix_type adsorb_over,adsorb_ham; /* overlap and hamiltonian matrix for adsorbate */
  eigenset_type adeigenset;
  real e_old,*save_e;
  wavepacket_type projection;
  real norm,normad;  /* For normalizing eigenvectors */
  real time;         /* Time during dynamics */
  BOOLEAN ad1,ad2;
  int num_adorbs;
  int istep;
  static char cubefile[100],wavefile[100],restartfile[100];
  char instring[200];
  FILE *f,*f_restart;
  
  fprintf(output_file,"\n# ************ELECTRON DYNAMICS*******************\n");
  
  /* 
   * Allocate memory for wavepacket  
   * */  
  ao.dim = num_orbs;
  ao.re = (real *)calloc(ao.dim,sizeof(real));
  ao.im = (real *)calloc(ao.dim,sizeof(real));
  save_ao.dim = num_orbs;
  save_ao.re = (real *)calloc(ao.dim,sizeof(real));
  save_ao.im = (real *)calloc(ao.dim,sizeof(real));
  mo.dim = num_orbs;
  mo.re = (real *)calloc(mo.dim,sizeof(real));
  mo.im = (real *)calloc(mo.dim,sizeof(real));
  save_mo.dim = num_orbs;
  save_mo.re = (real *)calloc(mo.dim,sizeof(real));
  save_mo.im = (real *)calloc(mo.dim,sizeof(real));

  /*-----------------Prepare Initial Wavepacket---------------------------*/

  if (cell->dim>1){  // There are periodic boundary conditions

    /* Copy K-space overlap to R-space overlap matrix */
    for (i=0;i<num_orbs;i++){
      for (j=i+1;j<num_orbs;j++){
        overlapR.mat[num_orbs*i+j] = overlapK.mat[num_orbs*i+j];
        overlapR.mat[num_orbs*j+i] = overlapK.mat[num_orbs*j+i];
	
	hamilR.mat[num_orbs*i+j] = hamilK.mat[num_orbs*i+j];
	hamilR.mat[num_orbs*j+i] = hamilK.mat[num_orbs*j+i];
      }
      overlapR.mat[num_orbs*i+i] = overlapK.mat[num_orbs*i+i];
    }

    for (i=0;i<cell->num_atoms;i++){
      //printf("atom[%d] = %d\n",i+1,orbital_lookup_table[i]);
    }
    //exit(1);
  }

  
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
   * If there are absorbing potentials
   * create them from the atom list and values
   * */
  details->absorbing_pot = (real *)calloc(num_orbs,sizeof(real));
  for (i=0;i<num_orbs;i++){
    details->absorbing_pot[i] = 0.0;
  }
  if (details->absorbing){
    for (i=0;i<details->absorbing_dim;i++){
      j = details->absorbing_orbs[i];
      details->absorbing_pot[j] = details->absorbing_value;
    }
  }

  /* 
   * Create orbital list for following time-dependent occupation
   *  */
  details->occ_orbs = (real *)calloc(num_orbs,sizeof(real));
  for (i=0;i<num_orbs;i++){
    details->occ_orbs[i] = 0.0;
  }
  if (details->occupation){
    for (i=0;i<details->occ_dim;i++){
      k = details->occ_atoms[i];
      if (k+1>cell->num_atoms-1){
	l = num_orbs - orbital_lookup_table[k];  
      } else{
        l = orbital_lookup_table[k+1]-orbital_lookup_table[k];
      }
      for (j=0;j<l;j++){
	details->occ_orbs[orbital_lookup_table[k]+j] = 1.0;
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

  fprintf(output_file,"\nAdsorbate Energies\n");
  for (i=0;i<num_adorbs;i++){
    fprintf(output_file,"%d:---> %f\n",i+1,adeigenset.val[i]);
  }
  
  /*
  *  Save energies  
  */
  save_e = (double *)calloc(num_orbs,sizeof(real));
  for (i=0;i<num_orbs;i++){
    save_e[i] = ENERGY(eigenset,i);
  }
  /*
   * Initalize Wavepacket *
   * */
  for (i=0;i<ao.dim;i++){
    ao.re[i]=0.0;
    ao.im[i]=0.0;
    save_ao.re[i]=0.0;
    save_ao.im[i]=0.0;
    mo.re[i]=0.0;
    mo.im[i]=0.0;
    save_mo.re[i]=0.0;
    save_mo.im[i]=0.0;
  }
  if (details->restart != 1){
    /* Check that initial packet superposition is normalized */
    norm = 0.0;
    for (j=0;j<details->initial_len;j++){
      norm += pow(details->initial_super[j],2.0);
    }
    norm = pow(norm,-0.5);
    for (j=0;j<details->initial_len;j++){
      details->initial_super[j] *= norm;
    }
    /* Loop over initial packet superposition 
     *      */
    for (j=0;j<details->initial_len;j++){
      k = details->initial_orbital[j]-1;
      for (i=0;i<num_adorbs;i++){
        //printf("adsorb orb %d\n",details->adsorbate_orbs[i]);
        ao.re[details->adsorbate_orbs[i]] += details->initial_super[j]*MATRIX(adsorb_ham,k,i);
        save_ao.re[details->adsorbate_orbs[i]] = ao.re[details->adsorbate_orbs[i]];
      }
    }
  }
  else{ /* Read initial state from restart file */
    fprintf(output_file,"Reading inital packet from restart file %s\n",details->restart_name);
    f_restart = fopen(details->restart_name,"r");
    fgets(instring,200,f_restart);
    sscanf(instring,"%d",&j);  /* Read number of elements in file*/
    if (j!=num_orbs){
      printf("Number of elements in restart file does not match system size. Exiting!\n");
      exit(1);
    }
    for (i=0;i<num_orbs;i++){
      fgets(instring,200,f_restart);
      sscanf(instring,"%lf %lf %lf %lf %lf",&(ao.re[i]),&(ao.im[i]),&(mo.re[i]),&(mo.im[i]),&e_old);
      //ENERGY(eigenset,i) = 0.5*(ENERGY(eigenset,i)+e_old);
    }
    fclose(f_restart);

    /* Calculate Projection of molecular orbitals on wave function
     *  <Phi_k(t)|Psi(t)> */
    projection.dim = num_orbs;
    projection.re = (real *)calloc(num_orbs,sizeof(real));
    projection.im = (real *)calloc(num_orbs,sizeof(real));
    project(&projection,ao,eigenset,overlapR);
    //for (i=0;i<num_orbs;i++){
    //  printf("%3d %22.15f %22.15f\n",i+1,projection.re[i],projection.im[i]);
    //}

    /* Copy Projection to mo's */
    for (i=0;i<num_orbs;i++){
      mo.re[i] = projection.re[i];
      mo.im[i] = projection.im[i];
    }
    /* Free projection */
    free(projection.re);
    free(projection.im);

    /* Transform MO's -> AO's by multipling with eigenvectors */
    fprintf(output_file,"\nComparing old ao's with new ao's\n");
    for (i=0;i<num_orbs;i++){
      fprintf(output_file,"%12.6f %12.6f",ao.re[i],ao.im[i]);
      ao.re[i] = 0.0;
      ao.im[i] = 0.0;
      for (j=0;j<num_orbs;j++){
        ao.re[i] += TRANSFORM(eigenset,j,i)*mo.re[j];
        ao.im[i] += TRANSFORM(eigenset,j,i)*mo.im[j];
      }
      fprintf(output_file,"%12.6f %12.6f\n",ao.re[i],ao.im[i]);
    }        

    /* Save AOs for calculation of survival probability */
    for (i=0;i<num_orbs;i++){
      save_ao.re[i] = ao.re[i];
      save_ao.im[i] = ao.im[i];
    }
  }
  norm = normalization(ao,overlapR);
  printf("Is the initial wavepacket normalized? %f\n",norm);
  
  /* 
   * Transform Initial wavepacet from ao basis -> mo basis
   * */
  /* First Multiply by Overlap */
  for (i=0;i<num_orbs;i++){
    mo.re[i] = 0.0;
    for (j=0;j<num_orbs;j++){
      mo.re[i] += HERMETIAN_R(overlapR,i,j)*ao.re[j];
    }
    save_mo.re[i] = mo.re[i]; /* Save result for Second Multiply */
  }        
  /* Then Multiply by Eigenvectors */
  for (i=0;i<num_orbs;i++){
    mo.re[i] = 0.0;
    for (j=0;j<num_orbs;j++){
      mo.re[i] += TRANSFORM(eigenset,i,j)*save_mo.re[j];
    }
  }
  for (i=0;i<num_orbs;i++){
    save_mo.re[i] = mo.re[i]; /* Now save Result */
  }
     

  fprintf(output_file,"\nInitial Wavepacket\n AO Basis    MO Basis\n");
  for (i=0;i<num_orbs;i++){
    fprintf(output_file," %8.5f  %8.5f  Phi(%4d)  %8.5f\n",ao.re[i],mo.re[i],i+1,ENERGY(eigenset,i));
  }
  fprintf(output_file,"\nAbsorbing Potential\n");
  for (i=0;i<num_orbs;i++){
    fprintf(output_file," %8.5f  E(%4d)\n",details->absorbing_pot[i],i+1);
  }

  /* 
   * The system is now ready for dynamics!!!
   * */

  /*--------------Electron Dynamics----------------------------------*/
  f = fopen(file_name,"w");
  fprintf(f,"# Electron Dynamics\n# Time  Survival\n");
  
  istep = 0;
  time = 0.0;

  /* prepare wave filename */
  if (details->cube==1){
    sprintf(wavefile,"%s.wave",file_name);
    printf("The wave file is named %s\n",&wavefile);

    /* call function to write header for wavefile */
    wave_header(wavefile,cell,num_orbs,orbital_lookup_table);
  
  }
  /* write initial wavepacket to wave file */
  if (details->cube==1){
    wave_packet(wavefile,ao);
  }

  /* Prepare restart file name */
  sprintf(restartfile,"%s.restart",file_name);
  
  /* Write Initial Wavepacket Data to file */
  fprintf(f,"%f   %f  %f  %f\n",time,survival(ao,overlapR,details),normalization(ao,overlapR),
		    occupation(ao,overlapR,details));

  write_restart(restartfile,ao,mo,save_e);

  printf("|t=%3.2f> = ",details->time_length);
  while (time < details->time_length){

    /*
     *  Propagate AO 
     *  */
    propagate(eigenset,details,&ao,mo,details->timestep);
    
    /* 
     * Transform time evolved Wavepacet from ao basis -> mo basis
     * */
    /* First Multiply by Overlap */
    for (i=0;i<num_orbs;i++){
      save_ao.re[i] = 0.0; // Use save_ao as working space
      save_ao.im[i] = 0.0;
      for (j=0;j<num_orbs;j++){
        save_ao.re[i] += HERMETIAN_R(overlapR,i,j)*ao.re[j];
        save_ao.im[i] += HERMETIAN_R(overlapR,i,j)*ao.im[j];
      }
    }        
    /* Then Multiply by Eigenvectors */
    for (i=0;i<num_orbs;i++){
      mo.re[i] = 0.0;
      mo.im[i] = 0.0;
      for (j=0;j<num_orbs;j++){
        mo.re[i] += TRANSFORM(eigenset,i,j)*save_ao.re[j];
        mo.im[i] += TRANSFORM(eigenset,i,j)*save_ao.im[j];
      }
    }        

    time += details->timestep;

    if (normalization(ao,overlapR)>1.1){ 
      fprintf(output_file,"\nAbsorbing Potential\n");
      for (i=0;i<num_orbs;i++){
         fprintf(output_file," %8.5f  E(%4d)\n",details->absorbing_pot[i],i+1);
      }
      exit(1);
    }


    /* Write Time-evolved Wavepacket Data to file */
    fprintf(f,"%f   %f  %f  %f\n",time,survival(ao,overlapR,details),normalization(ao,overlapR),
		    occupation(ao,overlapR,details));

    /* write wavepacket to wave file */
    if (details->cube==1){
      wave_packet(wavefile,ao);
    }
    
    /*write restart file */
    write_restart(restartfile,ao,mo,save_e);
    
    istep += 1;

    /* Time loop finished */
  } 
  printf("exp^(-iHt/hbar)*|t=%3.2f>\n",0.0); 
  
  /*-------------Done with Dynamics----------------------------------*/
  fclose(f);
  free(save_e);
}

/************************ propagate ************************************/

/* Propagation in Atomic Orbital Basis - Used with Absorbing Potentials 
 *
 * See Equation 1,2 in JACS vol 125, NO 26, 2003
 * 
 * */
void propagate(eigen,details,ao,mo,time)
  eigenset_type eigen;
  detail_type *details;
  wavepacket_type *ao,mo;
  real time;
{
  real hbar=0.65821189916; // eV*fs (electronVolt femtosecond)
  real save_psi;
  int i,j,k;

  for (i=0;i<ao->dim;i++){
    ao->re[i] = 0.0;
    ao->im[i] = 0.0;
    for (j=0;j<mo.dim;j++){
      ao->re[i] += mo.re[j]*TRANSFORM(eigen,j,i)*cos(ENERGY(eigen,j)*time*(1/hbar))+
	      mo.im[j]*TRANSFORM(eigen,j,i)*sin(ENERGY(eigen,j)*time*(1/hbar));
      ao->im[i] += -mo.re[j]*TRANSFORM(eigen,j,i)*sin(ENERGY(eigen,j)*time*(1/hbar))+
	      mo.im[j]*TRANSFORM(eigen,j,i)*cos(ENERGY(eigen,j)*time*(1/hbar));
    }
    ao->re[i] *= exp(-details->absorbing_pot[i]*time*(1/hbar));
    ao->im[i] *= exp(-details->absorbing_pot[i]*time*(1/hbar));
  }
  /* The wavepacket evolved! */
}


/********************** survival *********************************/

/* 
 * See Equation 6 in JACS vol 125, NO 26, 2003
 * */
real survival(ao,overlapR,details)
  wavepacket_type ao;
  hermetian_matrix_type overlapR;
  detail_type *details;
{
  real surv;
  int i,j,k;
  
  /* Calculate <t=0|t=time> */
  surv=0.0;
  for (i=0;i<ao.dim;i++){
    for (k=0;k<details->adsorbate_num_orbs;k++){
      j = details->adsorbate_orbs[k];
      surv += (ao.re[i]*ao.re[j]+ao.im[i]*ao.im[j])*HERMETIAN_R(overlapR,i,j);  /* Remember AO is not orthogonal */
    }
  }
  return surv;
}

/*********************** normalization ************************/
real normalization(ao,overlapR)
  wavepacket_type ao;
  hermetian_matrix_type overlapR;
{
  real norm,inorm;
  int i,j;

  /* Calculate <psi|psi> */
  norm = 0.0;
  for (i=0;i<ao.dim;i++){
    norm += (ao.re[i]*ao.re[i]+ao.im[i]*ao.im[i])*MATRIX(overlapR,i,i);
  }    
  norm = 0.0;
  inorm = 0.0;
  for (i=0;i<ao.dim;i++){
    for (j=0;j<ao.dim;j++){
      norm += (ao.re[i]*ao.re[j]+ao.im[i]*ao.im[j])*HERMETIAN_R(overlapR,i,j);
      inorm += (ao.re[i]*ao.im[j]-ao.im[i]*ao.re[j])*HERMETIAN_R(overlapR,i,j);
    }
  }
  return norm;
}

/*********************** occupation ********************************/
real occupation(ao,overlapR,details)
  wavepacket_type ao;
  hermetian_matrix_type overlapR;
  detail_type *details;
{
  real occ;
  int i,j,k;
  
  /* Calculate */
  occ = 0.0;
  for (i=0;i<ao.dim;i++){
    for (j=0;j<ao.dim;j++){
       occ += (ao.re[i]*ao.re[j]+ao.im[i]*ao.im[j])*HERMETIAN_R(overlapR,i,j)*details->occ_orbs[j];
    }
  }
  return occ; 
}


/********************* wave header **********************************/

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

/**************************** wave packet **********************************/

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
    fprintf(f,"%16.12lf  %16.12lf\n",ao.re[i],ao.im[i]);
  }
  fprintf(f,"\n");

  /* Close wavefile */
  fclose(f);
}

/**************************** write restart *********************************/

void write_restart(restartfile,ao,mo,energy)
  char *restartfile;
  wavepacket_type ao,mo;
  real *energy;
{
  int i;
  FILE *f;
  f = fopen(restartfile,"w");

  fprintf(f,"%d\n",ao.dim);
  for (i=0;i<ao.dim;i++){
    fprintf(f,"%20.15f %20.15f %20.15f %20.15f %20.15f\n",ao.re[i],ao.im[i],mo.re[i],mo.im[i],
              energy[i]);
  }
  fclose(f);
}


void project(projection,ao,eigenset,overlapR)
  wavepacket_type *projection,ao;
  eigenset_type eigenset;
  hermetian_matrix_type overlapR;
{
  int i,j,k;

  for (i=0;i<projection->dim;i++){
    projection->re[i] = 0.0;
    projection->im[i] = 0.0;
    for (j=0;j<ao.dim;j++){
      for (k=0;k<ao.dim;k++){
        /* real part */
        projection->re[i] += ao.re[j]*TRANSFORM(eigenset,i,k)*HERMETIAN_R(overlapR,j,k);
        /* imaginary part */
        projection->im[i] += ao.im[j]*TRANSFORM(eigenset,i,k)*HERMETIAN_R(overlapR,j,k);
      }
    }
  }
}
 
