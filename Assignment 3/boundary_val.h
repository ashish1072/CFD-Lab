#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

#include<math.h>
/**
 * The boundary values of the problem are set.
 */
void spec_boundary_val(int imax,int jmax,  double U_INFLOW,
  double V_INFLOW,double **U,double **V,int **flag);

void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V,
  int **flag
);

int B_O(int flag);

int B_W(int flag);

int B_N(int flag);

int B_S(int flag);

int B_NO(int flag);

int B_NW(int flag);

int B_SO(int flag);

int B_SW(int flag);

void set_temp_boundary(int imax, int jmax, double TI, double** T, int** Flag);

void set_dirichlet(int x1start, int x1end,int y1start,int y1end, double** T_set,double** T,int** Flag);
#endif
