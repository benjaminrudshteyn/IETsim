#include "cube.h"


double orbital(int n,int ml,double coeff1,double exp1,double coeff2,
		double exp2,double *r,double *dx,double *dy,double *dz);
double radial(int n,double coef1,double exp1,double coef2,
		double exp2,double *r);
int fact(int n);
double s_orb(void);
double p_x(double x,double y,double z);
double p_y(double x,double y,double z);
double p_z(double x,double y,double z);
double d_x2y2(double x,double y,double z);
double d_z2(double x,double y,double z);
double d_yz(double x,double y,double z);
double d_xz(double x,double y,double z);
double d_xy(double x,double y,double z);

double voxel(mol,basis,wave,x,y,z)
  molecule_type *mol;
  basis_type *basis;
  wavepacket_type *wave;
  double x,y,z;
{
  int i,j,k;	
  double au = 1.889726125; /* convert from Angstrom to Atomic Units */
  double rad,ang,vox_value;
  double dx,dy,dz,r;
  double dx2,dy2,dz2,r2;
  double a,areal,aimag,psi;
  
  areal = 0.0;
  aimag = 0.0;
	 
  /* Loop through wave */
  for (i=0;i<wave->num_orbs;i++){

    dx = (basis->x[i]-x)*au;
    dy = (basis->y[i]-y)*au;
    dz = (basis->z[i]-z)*au;   
    r = pow(pow(dx,2)+pow(dy,2)+pow(dz,2),0.5); 
	
    psi = orbital(basis->n[i],basis->ml[i],
	    basis->coeff1[i],basis->exp1[i],basis->coeff2[i],
 	    basis->exp2[i],&r,&dx,&dy,&dz);
    areal += wave->re[i]*psi;
    aimag += wave->im[i]*psi;
  }  
  if (density){
    vox_value = pow(areal,2.0)+pow(aimag,2.0);
  } else {
    vox_value = areal;
  }
  return vox_value;
}


double orbital(n,ml,coeff1,exp1,coeff2,exp2,r,dx,dy,dz)
  int n,ml;
  double coeff1,exp1,coeff2,exp2;
  double *r,*dx,*dy,*dz;
{
  double rad,ang;

  rad = radial(n,coeff1,exp1,coeff2,exp2,r); 
  switch (ml){
  case 1:
    ang = s_orb();
    break;
  case 2:
    ang = p_x(*dx,*dy,*dz);
    break;
  case 3:
    ang = p_y(*dx,*dy,*dz);
    break;
  case 4:
    ang = p_z(*dx,*dy,*dz);
    break;
  case 5:
    ang = d_x2y2(*dx,*dy,*dz);
    break;
  case 6:
    ang = d_z2(*dx,*dy,*dz);
    break;
  case 7:
    ang = d_xy(*dx,*dy,*dz);
    break;
  case 8:
    ang = d_xz(*dx,*dy,*dz);
    break;
  case 9:
    ang = d_yz(*dx,*dy,*dz);
    break;
  }
  return rad*ang;
}
  

double radial(int n,double coef1,double exp1,double coef2,
		double exp2,double *r){
  double norm1,norm2,value;

  norm1 = pow(2*exp1,2*n+1.0);
  norm2 = fact(2*n);

  value = pow(norm1/norm2,0.5)*pow(*r,n-1.0)*(coef1*exp(-exp1*(*r)));

  norm1 = pow(2*exp2,2*n+1.0);
  norm2 = fact(2*n);
	  
  value += pow(norm1/norm2,0.5)*pow(*r,n-1.0)*(coef2*exp(-exp2*(*r)));
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

double s_orb(void){
  return 0.28209479177237814;
}

double p_z(double x,double y,double z){
  double r;
  r = pow(x*x+y*y+z*z,0.5);
  if (r==0.0){
    return 0.0;
  } else{
    return 0.48860251190291992*(z/r);
  }
}


double p_y(double x,double y,double z){
  double r;
  r = pow(x*x+y*y+z*z,0.5);
  if (r==0.0){
    return 0.0;
  } else{
    return 0.48860251190291992*(y/r);
  }
}


double p_x(double x,double y,double z){
  double r;
  r = pow(x*x+y*y+z*z,0.5);
  if (r==0.0){
    return 0.0;
  } else{
    return 0.48860251190291992*(x/r);
  }
}

double d_x2y2(double x,double y,double z){
  double r2;
  r2 = x*x+y*y+z*z;
  if (r2==0.0){
    return 0.0;
  } else{
    return 1.0925484305920792*((x*x-y*y)/r2);
  }
}


double d_z2(double x,double y,double z){
  double r2;
  r2 = x*x+y*y+z*z;
  if (r2==0.0){
    return 0.0;
  } else{
    return 0.31539156525252005*(3*(z*z/r2)-1);
  }
}


double d_yz(double x,double y,double z){
  double r2;
  r2 = x*x+y*y+z*z;
  if (r2==0.0){
    return 0.0;
  } else{
    return 2.1850968611841584*((y*z)/r2);
  }
}


double d_xz(double x,double y,double z){
  double r2;
  r2 = x*x+y*y+z*z;
  if (r2==0.0){
    return 0.0;
  } else{
    return 2.1850968611841584*((x*z)/r2);
  }
}


double d_xy(double x,double y,double z){
  double r2;
  r2 = x*x+y*y+z*z;
  if (r2==0.0){
    return 0.0;
  } else{
    return 2.1850968611841584*((x*y)/r2);
  }
}
