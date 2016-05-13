#include "cube.h"



void parse_infile(char *filename,cube_type *cube);


int main(int argc,char *argv[]){
  molecule_type *mol;
  basis_type *basis;
  wavepacket_type *wave;
  cube_type *cube;
  int i,j,k;
  int n;
  int wave_index,make;
  char wavefile[100],infile[100];
  char line[230];
  FILE *f;

  n = sizeof(line);
  /*
   * Make space for needed types
   * */
  mol = (molecule_type *)calloc(1,sizeof(molecule_type));
  basis = (basis_type *)calloc(1,sizeof(basis_type));
  wave = (wavepacket_type *)calloc(1,sizeof(wavepacket_type));
  cube = (cube_type *)calloc(1,sizeof(cube_type));
  
  /* 
   * Get information about cube boundary
   * */
  parse_infile(argv[1],cube);


  /* 
   * Determine number of voxels and make space 
   * */
  cube->xvox = (int)((cube->xmax-cube->xmin)/cube->dr) + 1;
  cube->yvox = (int)((cube->ymax-cube->ymin)/cube->dr) + 1;
  cube->zvox = (int)((cube->zmax-cube->zmin)/cube->dr) + 1;
  cube->num_vox = cube->xvox * cube->yvox * cube->zvox; 
  cube->vox = (double *)calloc(cube->num_vox,sizeof(double));
  if (debug){
     printf("\nCheck voxel data\n #xvox = %d #yvox = %d #zvox = %d\n",
		     cube->xvox,cube->yvox,cube->zvox);
  }
  
  /* 
   * Start Reading Wavefile 
   * */
  f = fopen(argv[2],"r");
  mol->num_atoms = 0;    /* set to false */
  basis->num_basis = 0;
  wave->num_orbs = 0;
  wave_index = 0;
  make = 0;
  while(nextline(f,&line,n)){
    /* make line lower case */
    for (i=0;i<strlen(line);i++){
      line[i] = tolower(line[i]);
    }
    /*----------Check for Geometry---------*/
    if (strstr(line,"molecule")){
      nextline(f,&line,n);
      if (!mol->num_atoms){
        sscanf(line,"%d",&(mol->num_atoms));
	mol->atom_number = (int *)calloc(mol->num_atoms,sizeof(int));
	mol->x = (double *)calloc(mol->num_atoms,sizeof(double));
	mol->y = (double *)calloc(mol->num_atoms,sizeof(double));
	mol->z = (double *)calloc(mol->num_atoms,sizeof(double));
      }
      for (i=0;i<mol->num_atoms;i++){
        nextline(f,&line,n);
	sscanf(line,"%d %lf %lf %lf",&(mol->atom_number[i]),
			&(mol->x[i]),&(mol->y[i]),&(mol->z[i]));
      }	      
    }
    /*----------Check for Basis------------*/
    else if (strstr(line,"basis")){
      nextline(f,&line,n);
      if (!basis->num_basis){
        sscanf(line,"%d",&(basis->num_basis));
	basis->n = (int *)calloc(basis->num_basis,sizeof(int));
	basis->ml = (int *)calloc(basis->num_basis,sizeof(int));
	basis->coeff1 = (double *)calloc(basis->num_basis,sizeof(double));
	basis->coeff2 = (double *)calloc(basis->num_basis,sizeof(double));
	basis->exp1 = (double *)calloc(basis->num_basis,sizeof(double));
	basis->exp2 = (double *)calloc(basis->num_basis,sizeof(double));
	basis->x = (double *)calloc(basis->num_basis,sizeof(double));
	basis->y = (double *)calloc(basis->num_basis,sizeof(double));
	basis->z = (double *)calloc(basis->num_basis,sizeof(double));
      }	   
      for (i=0;i<basis->num_basis;i++){
        nextline(f,&line,n);
	sscanf(line,"%d %d %d %lf %lf %lf %lf",&j,&(basis->n[i]),&(basis->ml[i]),
			&(basis->coeff1[i]),&(basis->exp1[i]),&(basis->coeff2[i]),
			&(basis->exp2[i]));
	basis->x[i] = mol->x[j];
	basis->y[i] = mol->y[j];
	basis->z[i] = mol->z[j];	
      }	    

      /*-------Check for Wavepacket----------*/
    }	 

    else if (strstr(line,"wave")){
      nextline(f,&line,n);
      if (!wave->num_orbs){
        sscanf(line,"%d",&(wave->num_orbs));
	wave->re = (double *)calloc(wave->num_orbs,sizeof(double));
	wave->im = (double *)calloc(wave->num_orbs,sizeof(double));
      }		
      for (i=0;i<wave->num_orbs;i++){
        nextline(f,&line,n);
        sscanf(line,"%lf %lf",&(wave->re[i]),&(wave->im[i]));	
      }
      /* Check if this cube is to be created */
      make = 0;
      for (i=0;i<cube->totnum_make;i++){
	if (wave_index==cube->make_which[i]) make = 1;
      }	   
      /* Build the cube, this is the hard work */      
      if (make) {
        printf("Building cube %d\n",wave_index);
        build(mol,basis,wave,cube,wave_index);
      }
      /* Build is done */
      make = 0;
      wave_index += 1;   
    }	    

  }	  
  fclose(f);  

  return 0;
}

void parse_infile(char *filename,cube_type *cube){
  int i,j,k;
  int start,stop;
  int m,n,o,base;
  int done;
  int section,factor,constant;
  char line[230];
  char *tok;
  char dash;
  FILE *f;

  f = fopen(filename,"r");

  n = sizeof(line);
  while(nextline(f,&line,n)){

    /* Make line uppercase */
    for (i=0;i<strlen(line);i++){
      line[i] = tolower(line[i]);
    }

    /* 
     * Look for keywords 
     * */

    /*-------------Cube Region-------------*/
    if (strstr(line,"region")){
      nextline(f,&line,n);
      sscanf(line,"%lf",&(cube->dr)); /* grid size */      
      nextline(f,&line,n);  /* x range */
      sscanf(line,"%lf %lf",&(cube->xmin),&(cube->xmax));
      nextline(f,&line,n);  /* y range */
      sscanf(line,"%lf %lf",&(cube->ymin),&(cube->ymax));
      nextline(f,&line,n);  /* z range */
      sscanf(line,"%lf %lf",&(cube->zmin),&(cube->zmax));
      if (debug){
	printf("\nReading Cube region\n");
	printf("dr = %f\n",cube->dr);
	printf("xmin = %f xmax = %f\n",cube->xmin,cube->xmax);
	printf("ymin = %f ymax = %f\n",cube->ymin,cube->ymax);
	printf("zmin = %f zmax = %f\n",cube->zmin,cube->zmax);
      }
    }

    /*-------------Which Cubes to make-----*/
    else if(strstr(line,"make")){
      /* Get number of sections */
      nextline(f,&line,n);
      sscanf(line,"%d %d",&section,&(cube->totnum_make));

      /* make space for index of which cubes will be made */
      cube->make_which = (int *)calloc(cube->totnum_make,sizeof(int));
      cube->num_make = (int *)calloc(section,sizeof(int));
  
      for (m=1;m<=section;m++){
        /* Get number of cube files to be made */
        nextline(f,&line,n);
        sscanf(line,"%d",&(cube->num_make[m]));
        nextline(f,&line,n);
        sscanf(line,"%d %d",&factor,&constant);
  
        /* read index of which cubes to make, seperated by commas */
        nextline(f,&line,n);
        tok = strtok(line,",\n");
        done = 0;
        j = 0;
        base = 0;
        if (m>1){
          for (o=1;o<m;o++){
            base = base + cube->num_make[o];
          }
        }
        while (!done){
          if (tok[strlen(tok)-1]=='\\'){  /* continue on next line */
            nextline(f,&line,n);
            tok = strtok(line,",\n");
          } else if (j>0) {
            tok = strtok(0,",\n");  /* next token */
          }
          if (!tok || tok[0]==0){   /* done */
            done = 1;
          } else {
            /* Check for dash */
            if (strstr(tok,"-")){
              printf("Found dash!\n");
              sscanf(tok,"%d%c%d",&start,&dash,&stop);
              for (k=start;k<=stop;k++){
                cube->make_which[j+base] = k*factor+constant;
                j += 1;
              }
            } else {
              sscanf(tok,"%d",&(cube->make_which[j+base])); /* store index */
              j += 1;
            }
          }
        }
        if (debug){
          printf("\nReading which cubes to make\n");
          for (j=0;j<cube->num_make[m];j++){
            printf("will make cube %d\n",cube->make_which[j+base]);
          }
        }
      }
    }

    /*-------------What to name cube files--*/
    else if(strstr(line,"name")){
      nextline(f,&line,n);
      sscanf(line,"%s",&(cube->name));
      if (debug){
        printf("\nWill name cubes: %s\n",&(cube->name));
      }
    }
  }

  /*
   * Done with input file
   * */

  fclose(f);
}

int nextline(FILE *f,char *line,int n){
  char noline[]="\n";
  char comment[]="#";
  int i;
  /* returns next line skipping comments (#) and blank lines
   * return 0 if EOF is reached */
  
  /* Get next line from file and Check for end of file */
  if(fgets(line,n*sizeof(char),f)==NULL) return 0;
  
  /* Check for blank line or comment */
  if(strcmp(line,noline)==0){
    return nextline(f,line,n); // recursion
  } else if(line[0]==comment[0]){
    return nextline(f,line,n); // recursion
  } else {
    return 1;  // base case
  }
}
