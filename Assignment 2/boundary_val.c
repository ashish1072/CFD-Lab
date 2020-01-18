 #include "boundary_val.h"
#include<string.h>

void boundaryvalues(int imax, int jmax, double** U,double **V,int **Flag)
{
  int i;
  int j;
	
	for(i = 0 ; i <= imax +1 ; i++ )
		for(j = 0 ; j <= jmax +1; j++ )
		{
			if(Flag[i][j] > 16)                                        
				switch(Flag[i][j]&480)
				{
					case 32:
									V[i][j] = 0;							//Boundary values for cell with B_N flag
									U[i-1][j] = -U[i-1][j+1];
									U[i][j] = -U[i][j+1];
									break;
					case 64:
									V[i][j-1] = 0;							//Boundary values for cell with B_S flag
									U[i-1][j] = -U[i-1][j-1];
									U[i][j] = -U[i][j-1];
									break;
					case 128:
									U[i-1][j] = 0;							//Boundary values for cell with B_W flag
									V[i][j-1] = -V[i-1][j-1];
									V[i][j] = -V[i-1][j];
									break;
					case 256:
									U[i][j] = 0;						    //Boundary values for cell with B_E flag
									V[i][j-1] = -V[i+1][j-1];
									V[i][j] = -V[i+1][j];
									break;
					case 288:
									V[i][j] = 0;	
									U[i][j] = 0;							//Boundary values for cell with B_NE flag
									U[i-1][j] = -U[i-1][j+1];
									V[i][j-1] = -V[i+1][j-1];
									break;
					case 192:
									V[i][j-1] = 0;	
									U[i-1][j] = 0;			                //Boundary values for cell with B_SW flag                                                                                                                     																	//Boundary values for cell with B_SW flag
									U[i][j] = -U[i][j-1];
									V[i][j] = -V[i-1][j];
									break;
					case 320:
									V[i][j-1] = 0;	
									U[i][j] = 0;							//Boundary values for cell with B_SE flag
									U[i-1][j] = -U[i-1][j-1];
									V[i][j] = -V[i+1][j];
									break;
					case 160:
									V[i][j] = 0;	
									U[i-1][j] = 0;						    //Boundary values for cell with B_NW flag
									U[i][j] = -U[i][j+1];
									V[i][j-1] = -V[i-1][j-1];
									break;





				}
		}	

} 


void spec_boundary_val( int imax, int jmax, double** U,double **V,int **Flag)
{
	for(int i = 0; i<= imax + 1; ++i)
	{
		for(int j = 0; j<= jmax +1; ++j)
		{
		    if(Flag[i][j] & 16)               // inflow condition
			{ 
			    U[i][j]=1; 
				V[i][j] = 0;         // *****
		        V[i][j-1] = 0;
			}

		    if(Flag[i][j] & 8)                // outflow condition(Right Boundary)
		    {
			U[i-1][j] = U[i-2][j];
			V[i][j] = V[i-1][j];
		    }
	    }
	}
		
	

}
