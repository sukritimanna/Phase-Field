// this is a routine which will take care of the periodic stresses

void avg_sigma_star(double average_periodic_stress[3][3],double kron_delta[3][3],fftw_complex *het11,fftw_complex *het12,fftw_complex *het22,fftw_complex *fi,double lambda_mat,double mu_mat,double mu_prime_mat,double del_lambda,double del_mu,double del_mu_prime,int nodes_x,int nodes_y,double dx,double dy,int iter)
{
	//double block_sigma[3][3],sigma_current[3][3],sigma_right[3][3],sigma_up[3][3],sigma_diagonal[3][3];
	double lambda,mu,mu_prime,het[3][3];
	int i,j,m,n;//,count=0;
	double sigma_current[3][3];
	long current;//,neighbour_right,neighbour_up,neighbour_diagonal;
	//initialising the average periodic stress matrix to zero
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			average_periodic_stress[i][j]=0.0;
		}
	}
	for(n=0;n<nodes_y;n++)
	{
		for(m=0;m<nodes_x;m++)
		{	
			//count++;
			current=m+nodes_x*n;
			//neighbour_right=(m+1)+nodes_x*n;
			//neighbour_up=m+nodes_x*(n+1);
			//neighbour_diagonal=(m+1)+nodes_x*(n+1);
			
			
			//populating the sigma_current tensor
			//calculating the value of lambda, mu and mu_prime for this node
			lambda=lambda_mat+del_lambda*creal(fi[current]);
			mu=mu_mat+del_mu*creal(fi[current]);
			mu_prime=mu_prime_mat+del_mu_prime*creal(fi[current]);
			//intialising the het matrix
			for(i=0;i<3;i++)
			{
				for(j=0;j<3;j++)
				{
					if(i==0 && j==0)
						het[i][j]=het11[current];
					else if(i==1 && j==1)
						het[i][j]=het22[current];
					else if(i==0 && j==1 || i==1 && j==0)
						het[i][j]=het12[current];
					else
						het[i][j]=0.0;
				}
			}
			//computing the stresses for this current node
			str(sigma_current,het,kron_delta,lambda,mu,mu_prime);
			
			//populating the sigma_right tensor
			//calculating the value of lambda, mu and mu_prime for this node
			//lambda=lambda_mat+del_lambda*creal(fi[neighbour_right]);
			//mu=mu_mat+del_mu*creal(fi[neighbour_right]);
			//mu_prime=mu_prime_mat+del_mu_prime*creal(fi[neighbour_right]);
			//intialising the het matrix
			//for(i=0;i<3;i++)
			//{
				//for(j=0;j<3;j++)
				//{
					//if(i==0 && j==0)
						//het[i][j]=het11[neighbour_right];
					//else if(i==1 && j==1)
						//het[i][j]=het22[neighbour_right];
					//else if((i==0 && j==1) || (i==1 && j==0))
						//het[i][j]=het12[neighbour_right];
					//else
						//het[i][j]=0.0;
				//}
			//}
			//computing the stresses for this current node
			//str(sigma_right,het,kron_delta,lambda,mu,mu_prime);
				
			//populating the sigma_up tensor
			//calculating the value of lambda, mu and mu_prime for this node
			//lambda=lambda_mat+del_lambda*creal(fi[neighbour_up]);
			//mu=mu_mat+del_mu*creal(fi[neighbour_up]);
			//mu_prime=mu_prime_mat+del_mu_prime*creal(fi[neighbour_up]);
			//intialising the het matrix
			//for(i=0;i<3;i++)
			//{
				//for(j=0;j<3;j++)
				//{
					//if(i==0 && j==0)
						//het[i][j]=het11[neighbour_up];
					//else if(i==1 && j==1)
						//het[i][j]=het22[neighbour_up];
					//else if(i==0 && j==1 || i==1 && j==0)
						//het[i][j]=het12[neighbour_up];
					//else
						//het[i][j]=0.0;
				//}
			//}
			//computing the stresses for this current node
			//str(sigma_up,het,kron_delta,lambda,mu,mu_prime);

			//populating the sigma_diagonal tensor
			//calculating the value of lambda, mu and mu_prime for this node
			//lambda=lambda_mat+del_lambda*creal(fi[neighbour_diagonal]);
			//mu=mu_mat+del_mu*creal(fi[neighbour_diagonal]);
			//mu_prime=mu_prime_mat+del_mu_prime*creal(fi[neighbour_diagonal]);
			//intialising the het matrix
			//for(i=0;i<3;i++)
			//{
				//for(j=0;j<3;j++)
				//{
					//if(i==0 && j==0)
					//	het[i][j]=het11[neighbour_diagonal];
					//else if(i==1 && j==1)
						//het[i][j]=het22[neighbour_diagonal];
					//else if(i==0 && j==1 || i==1 && j==0)
						//het[i][j]=het12[neighbour_diagonal];
					//else
						//het[i][j]=0.0;
				//}
			//}
			//computing the stresses for this current node
			//str(sigma_diagonal,het,kron_delta,lambda,mu,mu_prime);
			

			for(i=0;i<3;i++)
			{
				for(j=0;j<3;j++)
				{
					//computing the components of the block matrix 
					//block_sigma[i][j]=(sigma_current[i][j]+sigma_right[i][j]+sigma_up[i][j]+sigma_diagonal[i][j])/4.0; 
					//adding this contribution to the average_periodic_stress tensor
					average_periodic_stress[i][j]+=sigma_current[i][j];
				}
			}
		}
	}
	//dividing the summed over periodic stresses by the total system size
	//if(count!=((nodes_x-1)*(nodes_y-1)))
	//{
		//printf("Looping incorrect in average_periodic_stress.c routine\n");
		//exit(1);
	//}
	
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			average_periodic_stress[i][j]=average_periodic_stress[i][j]/(nodes_x*nodes_y);
			
		}
	
	} 
}
