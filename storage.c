//this is a routine to create the storage of the omega inverts 

void store(double omega_in[3][3],double *omega_in_store,int m,int n,int nodes_x)
{
	int i,k;
	for(i=0;i<3;i++)
	{
		for(k=0;k<3;k++)
		{ 
			omega_in_store[k+(3*i)+(3*3*m)+(3*3*nodes_x*n)]=omega_in[i][k];
			//this expression shows an equivalence between the double index notation and its single index variant 
		}
	}
}
