//this is a code for creating the omega and omega inverse matrices at every location and then storing up the inverse matrices 
void omega_cr(double *omega_in_store,int nodes_x,int nodes_y,double *kx,double *ky,double *sq_mag_k_vectors,double lambda_mat,double mu_mat,double mu_prime_mat,double kron_delta[3][3])
{
	double omega[3][3],omega_in[3][3],k_vector[3],max_omega;
	int m,n,i,k;
	long array_index;
	//creating the storage of the omega-invert matrices which are the LHS factors in the equations we are going to solve to get the heterogeneous elasticity
	for(n=0;n<nodes_y;n++)
	{
		for(m=0;m<nodes_x;m++)
		{
			array_index=m+nodes_x*n;
			//initialising the k_vector array with the values from the kx and ky arrays corresponding to that grid spacing
			//we consistently correlated the dimension indices 1,2 and 3 with x,y and z as in definitions of strain
			k_vector[0]=kx[m];
			k_vector[1]=ky[n];
			k_vector[2]=0.0;
			//printf("The omega matrix\n");
			//at the (0,0) grid location the k1,k2 and k3 are all zero so the omega matrix formed is not invertible   
			if(m!=0 || n!=0)
			{
				max_omega=0.0;
				for(i=0;i<3;i++)
				{
					for(k=0;k<3;k++)
					{
						//this expression considers the negative sign
						omega[i][k]=-(lambda_mat*k_vector[i]*k_vector[k]+mu_mat*kron_delta[i][k]*sq_mag_k_vectors[array_index]+mu_mat*k_vector[i]*k_vector[k]+mu_prime_mat*(k_vector[0]*k_vector[0]*kron_delta[i][0]*kron_delta[0][k]+k_vector[1]*k_vector[1]*kron_delta[i][1]*kron_delta[1][k]+k_vector[2]*k_vector[2]*kron_delta[i][2]*kron_delta[2][k]));	
						//printf("%lf\t",omega[i][k]);
					}
					//printf("\n");
				}
				//getchar();
				//determining the maximum of all the elements in the omega matrix
				for(i=0;i<3;i++)
				{
					for(k=0;k<3;k++)
					{
						if(max_omega<fabs(omega[i][k]))
							max_omega=fabs(omega[i][k]);
					}
				}
				//scaling down the omega matrix before inversion				 
				for(i=0;i<3;i++)
				{
					for(k=0;k<3;k++)
					{
						omega[i][k]=omega[i][k]/max_omega;
					}
				}
				//sending the matrix for inversion with pre-scaling (see "matrix_inversion.c")
			
				matrix_invert(omega,omega_in);
				//rescaling the omega_in matrix
				for(i=0;i<3;i++)
				{
					for(k=0;k<3;k++)
					{
						omega_in[i][k]=omega_in[i][k]/max_omega;
					}
				} 
			}
			else
			{
				//as omega-in is not defined 
				//populating the omega invert matrix with zeroesin the event of location being (0,0)
				//its implication is that the constant terms will be zero in the fourier representation of the displacements, i.e., we are not considering rigid body motion??? 
				for(i=0;i<3;i++)
				{
					for(k=0;k<3;k++)
					{
						omega_in[i][k]=0.0;
					}
				}
			}
			//sending the omega-invert matrix for storage (see "storage.c") 	
			store(omega_in,omega_in_store,m,n,nodes_x);
		}
	}
}
