#include "bind.h"

real radial(int n,real coef1,real exp1,real coef2,real exp2,real x,real y,real z);
real s_orb(void);
real p_x(real x,real y,real z);
real p_y(real x,real y,real z);
real p_z(real x,real y,real z);
real d_x2y2(real x,real y,real z);
real d_z2(real x,real y,real z);
real d_xy(real x,real y,real z);
real d_xz(real x,real y,real z);
real d_yz(real x,real y,real z);
int fact(int n);

/*      write_cube(cubefile,ao,overlapR,cell,orbital_lookup_table); */
void write_cube(filename,ao,overlapR,cell,orbital_lookup_table,details)
  char *filename;
  wavepacket_type ao;
  hermetian_matrix_type overlapR;
  cell_type *cell;
  int *orbital_lookup_table;
  detail_type *details;
{
/*
  Cube format:

  comment line      (ignored at import)
  comment line      (ignored at import)
  N  vx vy vz       number of atoms, origin of the volumetric data
  M1 vx1 vy1 vx1    number of voxels along first axis,
                    followed by first axis vector
  M2 vx2 vy2 vz2    number of voxels along second axis, second axis vector
  M2 vx3 vy3 vz3    numer of voxels, third axis vector
  atom1 x y z       atom number, coordinates
  atom2 x y z
  ...
  atomN x y z
  volumetric data   usually six voxes per line
*/
  int i,j,k,l;
  int ix,iy,iz;
  int iatom,jatom;
  int iorb,jorb;
  int iiorb,jjorb;
  int iorb_num,jorb_num;
  int xsteps,ysteps,zsteps;
  int num_atoms;
  real voxel,rad_iorb,ang_iorb,rad_jorb,ang_jorb;
  real zero = 0.0;
  real x,y,z;
  real xmax,xmin,ymax,ymin,zmax,zmin;
  real dx,dy,dz;
  real dr;  // step size along each axis
  real au = 1.889726125; //convert from Angstrom to Atomic Units
  FILE *f;
  f = fopen(filename,"w");
  
  dr = details->cube_step;
  xmin = details->cube_region[0];
  xmax = details->cube_region[1];
  ymin = details->cube_region[2];
  ymax = details->cube_region[3];
  zmin = details->cube_region[4];
  zmax = details->cube_region[5];
  
  /* 
   * Write comment lines 
   * */
  fprintf(f,"ROB CUBE FILE\nLOOP: X - Y- Z\n");

  /* determine number of atoms in cube region */
  num_atoms = 0;
  for (i=0;i<cell->num_atoms;i++){
    if (cell->atoms[i].loc.x>xmin && cell->atoms[i].loc.x<xmax){
      if (cell->atoms[i].loc.y>ymin && cell->atoms[i].loc.y<ymax){
	if (cell->atoms[i].loc.z>zmin && cell->atoms[i].loc.z<zmax){
          num_atoms += 1;
	}
      }
    }
  } 
  /* 
   * Write number of atoms and origin 
   * */
  fprintf(f,"%d  %f  %f  %f\n",num_atoms,xmin,ymin,zmin);

  /*
   *  Write Voxel vectors 
   *  */
  // First Vector - find number of voxels along axis
  xsteps = (int)((xmax-xmin)/dr) + 1;
  fprintf(f,"%d  %f  %f  %f\n",xsteps,dr,zero,zero);
  // Second Vector
  ysteps = (int)((ymax-ymin)/dr) + 1;
  fprintf(f,"%d  %f  %f  %f\n",ysteps,zero,dr,zero);
  // Third Vector
  zsteps = (int)((zmax-zmin)/dr) + 1;
  fprintf(f,"%d  %f  %f  %f\n",zsteps,zero,zero,dr);
  
  /* 
   * Write Molecule data to cube file
   * */
  /* write atoms and coordinates */
  for (i=0;i<cell->num_atoms;i++){
    if (cell->atoms[i].loc.x>xmin && cell->atoms[i].loc.x<xmax){
      if (cell->atoms[i].loc.y>ymin && cell->atoms[i].loc.y<ymax){
	if (cell->atoms[i].loc.z>zmin && cell->atoms[i].loc.z<zmax){
          fprintf(f,"%d  %f  %f  %f  %f\n",cell->atoms[i].at_number,
            zero,cell->atoms[i].loc.x,cell->atoms[i].loc.y,cell->atoms[i].loc.z);
	}
      }
    }
  }

  /* HARD PART!!
   *
   * Create Voxel data - loop over cartesian coords -> loop over orbitals
   *
   * This is incredibly nested and a more efficient algorithm should be found
   *  */
  
  /* Loop through volume */ 
  for (ix=0;ix<xsteps;ix++){
    x = (xmin+ix*dr);
    for (iy=0;iy<ysteps;iy++){
      y = (ymin+iy*dr);
      for (iz=0;iz<zsteps;iz++){
        z = (zmin+iz*dr);
        voxel = 0.0;

	/* Loop over atoms */
        for (iatom=0;iatom<cell->num_atoms;iatom++){
	  
          /* determine number of orbitals for this atom */
          if (iatom+1==cell->num_atoms){
	    iorb_num = overlapR.dim-orbital_lookup_table[iatom];
	  } else {
            iorb_num = orbital_lookup_table[iatom+1]-orbital_lookup_table[iatom];
	  }

	  /* determine distance of voxel from atom center */
          dx = (x-cell->atoms[iatom].loc.x)*au;
	  dy = (y-cell->atoms[iatom].loc.y)*au;
	  dz = (z-cell->atoms[iatom].loc.z)*au;

	  /*loop over orbitals for given atom (iatom) */
	  for (iorb=0;iorb<iorb_num;iorb++){
	    /* determine radial and angular contribution from this orbital */
	    switch (iorb){
            case 0:
	      rad_iorb = radial(cell->atoms[iatom].ns,1.0,cell->atoms[iatom].exp_s,0.0,0.0,dx,dy,dz);
	      ang_iorb = s_orb();
	      break;
            case 1:
	      rad_iorb = radial(cell->atoms[iatom].np,1.0,cell->atoms[iatom].exp_p,0.0,0.0,dx,dy,dz);
	      ang_iorb = p_x(dx,dy,dz);
	      break;
	    case 2:
	      rad_iorb = radial(cell->atoms[iatom].np,1.0,cell->atoms[iatom].exp_p,0.0,0.0,dx,dy,dz);
	      ang_iorb = p_y(dx,dy,dz);
	      break;
	    case 3:
	      rad_iorb = radial(cell->atoms[iatom].np,1.0,cell->atoms[iatom].exp_p,0.0,0.0,dx,dy,dz);
	      ang_iorb = p_z(dx,dy,dz);
	      break;
	    case 4:
	      rad_iorb = radial(cell->atoms[iatom].nd,cell->atoms[iatom].coeff_d1,cell->atoms[iatom].exp_d,
			      cell->atoms[iatom].coeff_d2,cell->atoms[iatom].exp_d2,dx,dy,dz);
	      ang_iorb = d_x2y2(dx,dy,dz);
	      break;
	    case 5:
	      rad_iorb = radial(cell->atoms[iatom].nd,cell->atoms[iatom].coeff_d1,cell->atoms[iatom].exp_d,
			      cell->atoms[iatom].coeff_d2,cell->atoms[iatom].exp_d2,dx,dy,dz);
	      ang_iorb = d_z2(dx,dy,dz);
	      break;
	    case 6:
	      rad_iorb = radial(cell->atoms[iatom].nd,cell->atoms[iatom].coeff_d1,cell->atoms[iatom].exp_d,
			      cell->atoms[iatom].coeff_d2,cell->atoms[iatom].exp_d2,dx,dy,dz);
	      ang_iorb = d_xy(dx,dy,dz);
	      break;
	    case 7:
	      rad_iorb = radial(cell->atoms[iatom].nd,cell->atoms[iatom].coeff_d1,cell->atoms[iatom].exp_d,
			      cell->atoms[iatom].coeff_d2,cell->atoms[iatom].exp_d2,dx,dy,dz);
	      ang_iorb = d_xz(dx,dy,dz);
	      break;
            case 8:
	      rad_iorb = radial(cell->atoms[iatom].nd,cell->atoms[iatom].coeff_d1,cell->atoms[iatom].exp_d,
			      cell->atoms[iatom].coeff_d2,cell->atoms[iatom].exp_d2,dx,dy,dz);
	      ang_iorb = d_yz(dx,dy,dz);
	      break;
	    }
	    /* loop over atoms again */
	    for (jatom=0;jatom<cell->num_atoms;jatom++){

              /* determine number of orbitals for this atom */
              if (jatom+1==cell->num_atoms){
     	        jorb_num = overlapR.dim-orbital_lookup_table[jatom];
	      } else {
                jorb_num = orbital_lookup_table[jatom+1]-orbital_lookup_table[jatom];
	      }

	      /* determine distance of voxel from atom center */
              dx = (x-cell->atoms[jatom].loc.x)*au;
	      dy = (y-cell->atoms[jatom].loc.y)*au;
	      dz = (z-cell->atoms[jatom].loc.z)*au;

              /* loop over orbitals for given atom (jatom) */
              for (jorb=0;jorb<jorb_num;jorb++){
	        switch (jorb){
                case 0:
	          rad_jorb = radial(cell->atoms[jatom].ns,1.0,cell->atoms[jatom].exp_s,0.0,0.0,dx,dy,dz);
	          ang_jorb = s_orb();
	          break;
		case 1:
	          rad_jorb = radial(cell->atoms[jatom].np,1.0,cell->atoms[jatom].exp_p,0.0,0.0,dx,dy,dz);
	          ang_jorb = p_x(dx,dy,dz);
	          break;
		case 2:
	          rad_jorb = radial(cell->atoms[jatom].np,1.0,cell->atoms[jatom].exp_p,0.0,0.0,dx,dy,dz);
	          ang_jorb = p_y(dx,dy,dz);
	          break;
		case 3:
	          rad_jorb = radial(cell->atoms[jatom].np,1.0,cell->atoms[jatom].exp_p,0.0,0.0,dx,dy,dz);
	          ang_jorb = p_z(dx,dy,dz);
	          break;
		case 4:
	          rad_jorb = radial(cell->atoms[jatom].nd,cell->atoms[jatom].coeff_d1,cell->atoms[jatom].exp_d,
			      cell->atoms[jatom].coeff_d2,cell->atoms[jatom].exp_d2,dx,dy,dz);
	          ang_jorb = d_x2y2(dx,dy,dz);
	          break;
		case 5:
	          rad_jorb = radial(cell->atoms[jatom].nd,cell->atoms[jatom].coeff_d1,cell->atoms[jatom].exp_d,
			      cell->atoms[jatom].coeff_d2,cell->atoms[jatom].exp_d2,dx,dy,dz);
	          ang_jorb = d_z2(dx,dy,dz);
	          break;
		case 6:
	          rad_jorb = radial(cell->atoms[jatom].nd,cell->atoms[jatom].coeff_d1,cell->atoms[jatom].exp_d,
			      cell->atoms[jatom].coeff_d2,cell->atoms[jatom].exp_d2,dx,dy,dz);
	          ang_jorb = d_xy(dx,dy,dz);
	          break;
		case 7:
	          rad_jorb = radial(cell->atoms[jatom].nd,cell->atoms[jatom].coeff_d1,cell->atoms[jatom].exp_d,
			      cell->atoms[jatom].coeff_d2,cell->atoms[jatom].exp_d2,dx,dy,dz);
	          ang_jorb = d_xz(dx,dy,dz);
	          break;
		case 8:
	          rad_jorb = radial(cell->atoms[jatom].nd,cell->atoms[jatom].coeff_d1,cell->atoms[jatom].exp_d,
			      cell->atoms[jatom].coeff_d2,cell->atoms[jatom].exp_d2,dx,dy,dz);
	          ang_jorb = d_yz(dx,dy,dz);
	          break;
	        }
	        /* 
	        * Now append add contribution from to voxel from orbital i and j 
	        * Prob(r) = Psi_i^*(r)Psi_j(r)=\sum_i \sum_j a_i^* a_j phi_i phi_j
	        * */
	        iiorb = orbital_lookup_table[iatom]+iorb;
	        jjorb = orbital_lookup_table[jatom]+jorb;
	        //voxel += (ao.re[iiorb]*ao.re[jjorb]+ao.im[iiorb]*ao.im[jjorb])*(rad_iorb*ang_iorb*rad_jorb*ang_jorb); 
	      }
	    }
	    iiorb = orbital_lookup_table[iatom]+iorb;
	    voxel += (ao.re[iiorb])*(rad_iorb*ang_iorb);  
	  }
        }
        fprintf(f,"%f  ",voxel);
        if (iz%6==5) fprintf(f,"\n");
      }
      fprintf(f,"\n");
    }
  }
  // The End
  fclose(f);
} 


real radial(int n,real coef1,real exp1,real coef2,real exp2,real x,real y,real z){
  real r,norm1,norm2,value;
  
  r = pow(x*x+y*y+z*z,0.5);

  norm1 = pow(2*exp1,2*n+1.0);
  norm2 = fact(2*n);

  value = pow(norm1/norm2,0.5)*pow(r,n-1.0)*(coef1*exp(-exp1*r));

  norm1 = pow(2*exp2,2*n+1.0);
  norm2 = fact(2*n);
	  
  value += pow(norm1/norm2,0.5)*pow(r,n-1.0)*(coef2*exp(-exp2*r));
  return value;
}

int fact(int n){
  switch (n){
  case 0:
    return 1;
    break;
  case 1:
    return 1;
    break;
  default:
    return n*fact(n-1);
    break;
  }
}

real s_orb(void){
  return pow(4*acos(-1.0),-0.5);
}

real p_z(real x,real y,real z){
  real r;
  r = pow(x*x+y*y+z*z,0.5);
  if (r==0.0){
    return 0.0;
  } else{
    return pow(3.0/(4*acos(-1.0)),0.5)*(z/r);
  }
}


real p_y(real x,real y,real z){
  real r;
  r = pow(x*x+y*y+z*z,0.5);
  if (r==0.0){
    return 0.0;
  } else{
    return pow(3.0/(4*acos(-1.0)),0.5)*(y/r);
  }
}


real p_x(real x,real y,real z){
  real r;
  r = pow(x*x+y*y+z*z,0.5);
  if (r==0.0){
    return 0.0;
  } else{
    return pow(3.0/(4*acos(-1.0)),0.5)*(x/r);
  }
}

real d_x2y2(real x,real y,real z){
  real r2;
  r2 = x*x+y*y+z*z;
  if (r2==0.0){
    return 0.0;
  } else{
    return pow(5.0/(16*acos(-1.0)),0.5)*((x*x-y*y)/r2);
  }
}


real d_z2(real x,real y,real z){
  real r2;
  r2 = x*x+y*y+z*z;
  if (r2==0.0){
    return 0.0;
  } else{
    return pow(5.0/(16*acos(-1.0)),0.5)*((z*z-r2)/r2);
  }
}


real d_yz(real x,real y,real z){
  real r2;
  r2 = x*x+y*y+z*z;
  if (r2==0.0){
    return 0.0;
  } else{
    return pow(5.0/(4*acos(-1.0)),0.5)*((y*z)/r2);
  }
}


real d_xz(real x,real y,real z){
  real r2;
  r2 = x*x+y*y+z*z;
  if (r2==0.0){
    return 0.0;
  } else{
    return pow(5.0/(4*acos(-1.0)),0.5)*((x*z)/r2);
  }
}


real d_xy(real x,real y,real z){
  real r2;
  r2 = x*x+y*y+z*z;
  if (r2==0.0){
    return 0.0;
  } else{
    return pow(5.0/(4*acos(-1.0)),0.5)*((x*y)/r2);
  }
}
