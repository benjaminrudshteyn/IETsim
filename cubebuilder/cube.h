/* Include file for the cubebuilder program */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define debug 1
#define density 1

/* MOLECULE Definition */
typedef struct {
  int num_atoms;
  int *atom_number;
  double *x;
  double *y;
  double *z;
} molecule_type;

/* Basis Definition */
typedef struct {
  int num_basis;
  int *n;
  int *ml;
  double *coeff1;
  double *coeff2;
  double *exp1;
  double *exp2;
  double *x;
  double *y;
  double *z;
} basis_type;

/* Wavepacket Definition */
typedef struct {
  int num_orbs;
  double *re;
  double *im;
} wavepacket_type;

/* Cube Definition */
typedef struct {
  char name[100]; /* name */

  /* region */
  double dr;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;

  /* voxel */
  int xvox;
  int yvox;
  int zvox;
  int num_vox;
  double *vox;

  /* make */
  int totnum_make;
  int *num_make;
  int *make_which;
} cube_type;
