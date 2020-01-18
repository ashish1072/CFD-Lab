#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V,
  int **Flag
);

void spec_boundary_val( int imax, int jmax, double** U,double **V,int **Flag);

#endif
