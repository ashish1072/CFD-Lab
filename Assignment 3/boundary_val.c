#include "boundary_val.h"
#include <stdio.h>

void spec_boundary_val(int imax,int jmax, double U_INFLOW, double V_INFLOW,double **U,double **V,int **flag)
{
	for(int i = 0; i<imax; ++i)
		for(int j = 0; j<jmax; ++j)
		    if(flag[i][j]&(1<<4))
			{ 

			    U[i][j]=U_INFLOW; V[i][j] = V_INFLOW; 
			    
			}
		
}
	
int B_O(int flag)
{
	return ( (flag & (1<<8)) && ~( flag & ((1<<5) | (1<<6)) ) );
}

int B_W(int flag)
{
	return ( (flag & (1<<7)) && ~( flag & ((1<<5) | (1<<6)) ) );
}

int B_N(int flag)
{
	return ( (flag & (1<<5)) && ~( flag & ((1<<7) | (1<<8)) ) );
}

int B_S(int flag) 
{
	return ( (flag & (1<<6)) && ~( flag & ((1<<7) | (1<<8)) ) );
}

int B_NO(int flag)
{
	return ( (flag & (1<<8)) &&  (flag & (1<<5)) );
}

int B_NW(int flag)
{
	return ( (flag & (1<<7)) &&  (flag & (1<<5)) );
}

int B_SO(int flag)
{
	return ( (flag & (1<<8)) &&  (flag & (1<<6)) ); 
}

int B_SW(int flag)
{
	return ( (flag & (1<<7)) &&  (flag & (1<<6)) );
}


void boundaryvalues(int imax,int jmax,double **U,double **V,int **flag)
{
	for(int i = 0; i<imax; ++i){
	    for(int j = 0; j<jmax; ++j){
			switch(flag[i][j]& ((1<<1)|(1<<2)|(1<<3)|(1<<9)) )
            {
            	case 1<<1://No slip conditions
		//					printf("It is no slip boundary flag[%d][%d] = %d\n",i,j,flag[i][j]);
            	            if ( B_O(flag[i][j]) )
                        	{
                            	U[i][j] = 0;
                            	V[i][j-1] = -V[i+1][j-1];		
                            	V[i][j] = -V[i+1][j];
                        	}
            
                        	if ( B_W(flag[i][j]) ) 
                        	{
                            	U[i-1][j] = 0;
                            	V[i][j-1] = -V[i-1][j-1];	
                            	V[i][j] = -V[i-1][j];
                        	}
                        
                        	if ( B_N(flag[i][j]) ) 
                        	{
                            	V[i][j] = 0;
                            	U[i-1][j] = -U[i-1][j+1];
                            	U[i][j] = -U[i][j+1];	
                        	}
                        
                        	if ( B_S(flag[i][j]) )
                        	{
                            	V[i][j-1] = 0;
                            	U[i-1][j] = -U[i-1][j-1];
                            	U[i][j] = -U[i][j-1];	
                        	}
                        
                        	if ( B_NO(flag[i][j])  )
                        	{
                            	U[i][j] = 0;
                            	U[i-1][j] = -U[i-1][j+1];
                            	V[i][j] = 0;
                            	V[i][j-1] = -V[i+1][j-1];
                        	}
                        
                        	if ( B_NW(flag[i][j])  )
                        	{
                            	U[i-1][j] = 0;
                            	U[i][j] = - U[i][j+1];
                            	V[i][j] = 0;
                            	V[i][j-1] = -V[i-1][j-1];
                        	} 
                        
                        	if ( B_SO(flag[i][j]) )
                        	{
                            	U[i][j]=0;
                            	U[i-1][j] = -U[i-1][j-1];
                            	V[i][j-1]=0;
                            	V[i][j] = -V[i+1][j];
                        	}
                        
                        	if ( B_SW(flag[i][j]) )
                        	{
                            	U[i-1][j] = 0;
                            	U[i][j] = -U[i][j-1];
                            	V[i][j-1] = 0;
                            	V[i][j] = -V[i-1][j];
                        	}
                        	break;
				case 1<<9://No slip conditions
		//				printf("It is coupling boundary flag[%d][%d] = %d\n",i,j,flag[i][j]);
            	            if ( B_O(flag[i][j]) )
                        	{
                            	U[i][j] = 0;
                            	V[i][j-1] = -V[i+1][j-1];		
                            	V[i][j] = -V[i+1][j];
                        	}
            
                        	if ( B_W(flag[i][j]) ) 
                        	{
                            	U[i-1][j] = 0;
                            	V[i][j-1] = -V[i-1][j-1];	
                            	V[i][j] = -V[i-1][j];
                        	}
                        
                        	if ( B_N(flag[i][j]) ) 
                        	{
                            	V[i][j] = 0;
                            	U[i-1][j] = -U[i-1][j+1];
                            	U[i][j] = -U[i][j+1];	
                        	}
                        
                        	if ( B_S(flag[i][j]) )
                        	{
                            	V[i][j-1] = 0;
                            	U[i-1][j] = -U[i-1][j-1];
                            	U[i][j] = -U[i][j-1];	
                        	}
                        
                        	if ( B_NO(flag[i][j])  )
                        	{
                            	U[i][j] = 0;
                            	U[i-1][j] = -U[i-1][j+1];
                            	V[i][j] = 0;
                            	V[i][j-1] = -V[i+1][j-1];
                        	}
                        
                        	if ( B_NW(flag[i][j])  )
                        	{
                            	U[i-1][j] = 0;
                            	U[i][j] = - U[i][j+1];
                            	V[i][j] = 0;
                            	V[i][j-1] = -V[i-1][j-1];
                        	} 
                        
                        	if ( B_SO(flag[i][j]) )
                        	{
                            	U[i][j]=0;
                            	U[i-1][j] = -U[i-1][j-1];
                            	V[i][j-1]=0;
                            	V[i][j] = -V[i+1][j];
                        	}
                        
                        	if ( B_SW(flag[i][j]) )
                        	{
                            	U[i-1][j] = 0;
                            	U[i][j] = -U[i][j-1];
                            	V[i][j-1] = 0;
                            	V[i][j] = -V[i-1][j];
                        	}
                        	break;
            
            	case 1<<2://Free Slip conditions
                        	if ( B_O(flag[i][j]) )
                        	{
                            	U[i][j] = 0;	
                            	V[i][j] = V[i+1][j];
                            	V[i][j-1] = V[i+1][j-1];
                        	}
                        
                        	if ( B_W(flag[i][j]) ) 
                        	{
                            	U[i-1][j] = 0;
                            	V[i][j-1] = V[i-1][j-1];	
                            	V[i][j] = V[i-1][j];
                        	}
                        
                        	if ( B_N(flag[i][j]) ) 
                        	{
                            	V[i][j] = 0;
                            	U[i-1][j] = U[i-1][j+1];
                            	U[i][j] = U[i][j+1];	
                        	}
                        
                        	if ( B_S(flag[i][j]) )
                        	{
                            	V[i][j-1] = 0;
                            	U[i-1][j] = U[i-1][j-1];
                            	U[i][j] = U[i][j-1];	
                        	}
                        
                        	if ( B_NO(flag[i][j])  )
                        	{
                            	U[i][j] = 0;
                            	U[i-1][j] = U[i-1][j+1];
                            	V[i][j] = 0;
                            	V[i][j-1] = V[i+1][j-1];
                        	}
                        
                        	if ( B_NW(flag[i][j])  )
                        	{
                            	U[i-1][j] = 0;
                            	U[i][j] = U[i][j+1];
                            	V[i][j] = 0;
                            	V[i][j-1] = V[i-1][j-1];
                        	} 
                        
                        	if ( B_SO(flag[i][j]) )
                        	{
                            	U[i][j]=0;
                            	U[i-1][j] = U[i-1][j-1];
                            	V[i][j-1]=0;
                            	V[i][j] = V[i+1][j];
                        	}
                        
                        	if ( B_SW(flag[i][j]) )
                        	{
                            	U[i-1][j] = 0;
                            	U[i][j] = U[i][j-1];
                            	V[i][j-1] = 0;
                            	V[i][j] = V[i-1][j];
                        	}	
                        	break;
            
            	case 1<<3:                                                 //Outflow from left 
              //          	printf("It is outflow boundary flag[%d][%d] = %d\n",i,j,flag[i][j]);
							if ( i==imax-1 && j>0){
							U[i-1][j] = U[i-2][j];
                        	V[i][j] = V[i-1][j];
							}
							else if (i==0 && j>0)
							{
							U[i][j] = U[i+1][j];
                        	V[i][j] = V[i+1][j];
				
							}
							else if (j==0 && i>0)
							{
							U[i][j] = U[i][j+1];
                        	V[i][j] = V[i][j+1];
                      
							}
							else if (j==jmax-1 && i>0)
							{
							U[i][j] = U[i][j-1];
                        	V[i][j-1] = V[i][j-2];
							}
                        	break;
            }
			}
			}
}

void set_temp_boundary(int imax, int jmax, double TI, double** T, int** Flag){
for (int i=0;i<imax;i++){
	for(int j=0;j<jmax;j++){
	//	
		if (((Flag[i][j]&(1<<0)))||((Flag[i][j]>>9))){  // coupling boundary or fluid
			continue;
		}
		else if(((Flag[i][j]&2)==2) || ((Flag[i][j]&4)==4)) // free slip or no slip
		{//printf("Flag[%d][%d]=%d\n",i,j,Flag[i][j]);
//			printf("setting the boundary for the flag[%d][%d] =%d\n",i,j,Flag[i][j]);
			switch (Flag[i][j]>>5) 
			{
			case 1: //north
//				printf("North\n");
				T[i][j] = T[i][j+1];
				break;
			case 2: //south
//				printf("South\n");
				T[i][j]=T[i][j-1];
				break;
			case 4: //west
//				printf("W\n");
				T[i][j]=T[i-1][j];
				break;
			case 8: // east
//				printf("E\n");
				T[i][j]=T[i+1][j];
				break;
			case 9: //North east
//				printf("N E\n");
				T[i][j] = (T[i][j+1] + T[i+1][j])/2;
				break;
			case 5: // north west
//				printf("N W\n");
				T[i][j] = (T[i][j+1] + T[i-1][j])/2;
				break;
			case 10: // South East
//				printf("S E\n");
				T[i][j] = (T[i][j-1] + T[i+1][j])/2;
				break;
			case 6: // South west
//				printf("S W\n");
				T[i][j] = (T[i][j-1] + T[i-1][j])/2;
				break;
			}
		}
		else if(Flag[i][j]&(1<<3)){ // outflow
	//		printf("Outflow boundary for the flag[%d][%d] =%d\n",i,j,Flag[i][j]);
//			printf("It is outflow boundary flag[%d][%d] = %d\n",i,j,Flag[i][j]);
							if ( i==imax-1 && j>0){
							T[i][j] = T[i-1][j];
							}
							else if (i==0 && j>0)
							{
							T[i][j] = T[i+1][j];
							}
							else if (j==0 && i>0)
							{
							T[i][j] = T[i][j+1];
							}
							else if (j==jmax-1 && i>0)
							{
							T[i][j] = T[i][j-1];
							}
			
		}
		else if(Flag[i][j]&(1<<4) ) //inflow
		{
			T[i][j]=TI;
		}
		
		
	}
}
}
void set_dirichlet(int x1start, int x1end,int y1start,int y1end, double** T_set,double** T,int** Flag){
	for (int i=x1start;i<=x1end;++i){
		for (int j=y1start;j<=y1end;++j){
		//	printf("i,j = %d,%d\n",i,j);
			//printf("Dirichlet's Flag[%d][%d]=%d\n",i,j,Flag[i][j]);
			if (Flag[i][j]>>9==1){ // no slip or free slip
		//		printf("Dirichlet's Flag[%d][%d]=%d\n",i,j,Flag[i][j]);
				switch (Flag[i][j]>>5)
				{
				case 17: // north
		//			printf("North case T_set[%d][%d] = %f\n",i,j,T_set[i][j]);
					T[i][j]=2*T_set[i][j] - T[i][j+1];
					break;
				case 18: //south
		//			printf("South case T_set[%d][%d] = %f\n",i,j,T_set[i][j]);
					T[i][j]=2*T_set[i][j] - T[i][j-1];
					break;
				case 20: //west
		//		printf("West case T_set[%d][%d] = %f\n",i,j,T_set[i][j]);
					T[i][j]=2*T_set[i][j] - T[i-1][j];
					break;
				case 24: // east
		//			printf("East case T_set[%d][%d] = %f\n",i,j,T_set[i][j]);
					T[i][j]=2*T_set[i][j]-T[i+1][j];
					break;
				case 25: //North east
		//		printf("NE case T_set[%d][%d] = %f\n",i,j,T_set[i][j]);
					T[i][j]=2*T_set[i][j] - T[i][j+1];
					T[i][j]=2*T_set[i][j]-T[i+1][j];
					break;
				case 21: // north west
		//		printf("NW case T_set[%d][%d] = %f\n",i,j,T_set[i][j]);
					T[i][j]=2*T_set[i][j] - T[i][j+1];
					T[i][j]=2*T_set[i][j] - T[i-1][j];
					break;
				case 26: // South East
		//		printf("SE case T_set[%d][%d] = %f\n",i,j,T_set[i][j]);
					T[i][j]=2*T_set[i][j] - T[i][j-1];
					T[i][j]=2*T_set[i][j]-T[i+1][j];
					break;
				case 22: // South west
		//		printf("SW T_set[%d][%d] = %f\n",i,j,T_set[i][j]);
					T[i][j]=2*T_set[i][j] - T[i][j-1];
					T[i][j]=2*T_set[i][j] - T[i-1][j];
					break;
								

				}
			}
		}
	}
}


