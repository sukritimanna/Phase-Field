//this is a routine for computing the average eigen stress in the domain
void avg_sigma_not(double average_eigen_stress[3][3],double sigma_c[3][3],double sigma_eta[3][3],double del_sigma_c[3][3],double del_sigma_eta[3][3],fftw_complex *theta_c,fftw_complex *theta_eta,fftw_complex *fi_times_theta_c,fftw_complex *fi_times_theta_eta,int nodes_x,int nodes_y,double dx,double dy)
{
	//double block_sigma[3][3],sigma_current[3][3],sigma_right[3][3],sigma_up[3][3],sigma_diagonal[3][3];
	int m,n,i,j;//,count=0;
	long current;//,neighbour_right,neighbour_up,neighbour_diagonal;
	double sigma_current[3][3];		
	//summing over the domain
	//initialising the average eigen stress matrix to zero
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			average_eigen_stress[i][j]=0.0;
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
			//populating the different sigma matrices
			for(i=0;i<3;i++)
			{
				for(j=0;j<3;j++)
				{
					//populating the sigma_current tensor
					sigma_current[i][j]=sigma_c[i][j]*creal(theta_c[current])+sigma_eta[i][j]*creal(theta_eta[current])+del_sigma_c[i][j]*creal(fi_times_theta_c[current])+del_sigma_eta[i][j]*creal(fi_times_theta_eta[current]);
					//sigma_current[i][j]=sigma_c[i][j]*creal(theta_c[current])+del_sigma_c[i][j]*creal(fi_times_theta_c[current]);	 
					//populating the sigma_right tensor
					//sigma_right[i][j]=sigma_c[i][j]*creal(theta_c[neighbour_right])+sigma_eta[i][j]*creal(theta_eta[neighbour_right])+del_sigma_c[i][j]*creal(fi_times_theta_c[neighbour_right])+del_sigma_eta[i][j]*creal(fi_times_theta_eta[neighbour_right]);	 
					//populating the sigma_up tensor
					//sigma_up[i][j]=sigma_c[i][j]*creal(theta_c[neighbour_up])+sigma_eta[i][j]*creal(theta_eta[neighbour_up])+del_sigma_c[i][j]*creal(fi_times_theta_c[neighbour_up])+del_sigma_eta[i][j]*creal(fi_times_theta_eta[neighbour_up]);	 
					//populating the sigma diagonal tensor
					//sigma_diagonal[i][j]=sigma_c[i][j]*creal(theta_c[neighbour_diagonal])+sigma_eta[i][j]*creal(theta_eta[neighbour_diagonal])+del_sigma_c[i][j]*creal(fi_times_theta_c[neighbour_diagonal])+del_sigma_eta[i][j]*creal(fi_times_theta_eta[neighbour_diagonal]);	
					//creating up the block_sigma tensor
					//block_sigma[i][j]=(sigma_current[i][j]+sigma_right[i][j]+sigma_up[i][j]+sigma_diagonal[i][j])/4.0;
					//adding this contribution to the average_eigen_stress tensor
					average_eigen_stress[i][j]+=sigma_current[i][j]; 
				}
			} 
		//closing the space loops
		}
	}
	//dividing the summed over eigen stresses by the total system size
	//if(count!=((nodes_x-1)*(nodes_y-1)))
	//{
	//	printf("Looping incorrect in average_eigen_stress.c routine\n");
	//	exit(1);
	//}
	
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			average_eigen_stress[i][j]=average_eigen_stress[i][j]/(nodes_x*nodes_y);
			
		}
		
	}
			
}
