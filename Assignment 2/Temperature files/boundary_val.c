#include<stdio.h>

void boundaryvalues(int imax, int jmax, double **U, double **V, double **T, double T_h, double T_c)
{
  int i;
  int j;

  // Boundary Conditions for U

  for ( j = 1 ; j <= jmax ; j++ ){
		U[0][j] = 0.0;                // left Boundary
		U[imax][j] = 0.0;             // Right Boundary
	}

	for( i = 1 ; i <= imax ; i++ ){
	    U[i][0] = -U[i][1];           // Bottom Boundary
		U[i][jmax+1] = - U[i][jmax];  // Top Boundary
	}

  // Boundary Conditions for V

	for ( j = 1 ; j <= jmax; j++ ){
	    V[0][j] = -V[1][j];           // left Boundary
		V[imax+1][j] = -V[imax][j];   // Right Boundary
	}
	
	for( i = 1 ; i <= imax ; i++ ){
		V[i][0] = 0.0;                // Bottom Boundary
		V[i][jmax] = 0.0;             // Top Boundary
	}

  // Boundary Conditions for T
  
    for ( j = 1 ; j <= jmax; j++ )
	{
	    T[0][j] = 2*T_h - T[1][j];           // left Boundary
		T[imax+1][j] = 2*T_c - T[imax][j];   // Right Boundary
	}
	
	for( i = 1 ; i <= imax ; i++ )
	{
		T[i][0] = T[i][1];                   // Bottom Boundary
		T[i][jmax+1] = T[i][jmax];           // Top Boundary
	}
}

/*  
    Boundary conditions for (d) example




    Boundary conditions for (e) example
    
	for( j = 1 ; j <= jmax ; j++ )
	{
		T[0][j] = T[1][j];                   //  Left Boundary
		T[imax+1][j] = T[imax][j];           //  Right Boundary
	}

    for ( i = 1 ; i <= imax; i++ )
	{
	    T[i][0] = 2*T_h - T[i][1];           // Bottom Boundary
		T[i][jmax+1] = 2*T_c - T[i][jmax];   // Top Boundary
	}

     */
     
    
  



