#include "cube.h"

double voxel(molecule_type *mol,basis_type *basis,wavepacket_type *wave,double x,double y,double z);

void build(mol,basis,wave,cube,cube_number)
  molecule_type *mol;
  basis_type *basis;
  wavepacket_type *wave;
  cube_type *cube;
  int cube_number;
{
  int i,j,k,l;
  int in_region;
  double x,y,z;
  double au = 1.889726125; /* convert from Angstrom to Atomic Units */
  char filename[100];
  FILE *f;

  if (cube_number<10){
    sprintf(filename,"%s_000%d.cube",cube->name,cube_number);
  }
  else if (cube_number<100){
    sprintf(filename,"%s_00%d.cube",cube->name,cube_number);
  }
  else if (cube_number<1000){  
    sprintf(filename,"%s_0%d.cube",cube->name,cube_number);
  }
  else if (cube_number<10000){
    sprintf(filename,"%s_%d.cube",cube->name,cube_number);
  }

  /*
   * Get voxel data
   */
  l = 0;
  for (i=0;i<cube->xvox;i++){
    for (j=0;j<cube->yvox;j++){
      for (k=0;k<cube->zvox;k++){
        x = cube->dr*i+cube->xmin;
	y = cube->dr*j+cube->ymin;
	z = cube->dr*k+cube->zmin;
	cube->vox[l] = voxel(mol,basis,wave,x,y,z);
	l += 1;
      }
    }
  }
  

  /*
   *   Cube format:
   *
   *   comment line      (ignored at import)
   *   comment line      (ignored at import)
   *   N  vx vy vz       number of atoms, origin of the volumetric data
   *   M1 vx1 vy1 vx1    number of voxels along first axis,
   *                     followed by first axis vector
   *   M2 vx2 vy2 vz2    number of voxels along second axis, second axis vector
   *   M2 vx3 vy3 vz3    numer of voxels, third axis vector
   *   atom1 x y z       atom number, coordinates
   *   atom2 x y z
   *    ...
   *   atomN x y z
   *   volumetric data   usually six voxes per line
   *                                             */
  
  f = fopen(filename,"w");
  /* 
   * Write comment lines
   *  */
  fprintf(f,"ROB CUBE FILE\nLOOP: X - Y- Z\n");

  /* determine number of atoms in cube region */
  in_region = 0;
  for (i=0;i<mol->num_atoms;i++){
    if (mol->x[i]>cube->xmin && mol->x[i]<cube->xmax && mol->y[i]>cube->ymin &&
		    mol->y[i]<cube->ymax && mol->z[i]>cube->zmin &&
		    mol->z[i]<cube->zmax) in_region += 1;
  }	  

  /* write number of atoms and origin */
  fprintf(f,"%d  %f  %f  %f\n",in_region,cube->xmin*au,cube->ymin*au,cube->zmin*au);
  
  /* write vectors */
  fprintf(f,"%d  %f  %f  %f\n",cube->xvox,cube->dr*au,0.0,0.0);
  fprintf(f,"%d  %f  %f  %f\n",cube->yvox,0.0,cube->dr*au,0.0);
  fprintf(f,"%d  %f  %f  %f\n",cube->zvox,0.0,0.0,cube->dr*au);
 
  /* write atoms and coords */
  for (i=0;i<mol->num_atoms;i++){
    if (mol->x[i]>cube->xmin && mol->x[i]<cube->xmax && mol->y[i]>cube->ymin &&
		    mol->y[i]<cube->ymax && mol->z[i]>cube->zmin &&
		    mol->z[i]<cube->zmax) {
      fprintf(f,"%d  %f  %f  %f  %f\n",mol->atom_number[i],0.0,
		      mol->x[i]*au,mol->y[i]*au,mol->z[i]*au);
    }		    
  }	  

  /* write voxel data */
  j = 0;
  for (i=0;i<cube->num_vox;i++){
     fprintf(f,"%f  ",cube->vox[i]);
     j += 1;
     if (j==6){
       fprintf(f,"\n");
       j = 0;
     }       
  }
	  
  
  fclose(f);
  
}
