// this is a routine for computing the elastic stresses
void elas_chem_pot(fftw_complex *het11,fftw_complex *het22,fftw_complex *het12,fftw_complex *fi,fftw_complex *theta_c,fftw_complex *theta_eta,double hom_strain[3][3],double kron_delta[3][3],double epsilon_c[3][3],double epsilon_eta[3][3],double epsilon_c_mag,double epsilon_eta_mag,double lambda_mat,double mu_mat,double mu_prime_mat,double del_lambda,double del_mu,double del_mu_prime, int nodes_x, int nodes_y,double A_z,double delta,double radius,double dx,fftw_complex *mu_el_comp,fftw_complex *mu_el_eta,double *el_str_11,double *el_str_22,double *el_str_12,double sigma_c[3][3],double sigma_eta[3][3],double del_sigma_c[3][3],double del_sigma_eta[3][3],double av_comp_alpha,double av_comp_beta,double av_eta, int t_steps)
{
	//char fname[300];
	//FILE *fp;
	int m,n,i,j;
	double het[3][3],epsilon_elastic[3][3],sigma_elastic[3][3],del_sigma_el[3][3];
	double lambda,mu,mu_prime;
	long array_index;

//opening the files for writing the elastic stress data (along the centre of the domain along the x-direction
		
	//sprintf(fname,"../DATA/elastic_stress_profile_ep_c_%.4lfep_eta_%.4lf_A_z_%.1lf_delta_%.1lf.dat",epsilon_c_mag,epsilon_eta_mag,A_z,delta);
	//if((fp=fopen(fname,"w"))==NULL)
	//{
	//	printf("Initial profile can't be written\n");
	//	exit(1);
	//}
	
	//calculating the elastic stresses and writing the files
	for(n=0;n<nodes_y;n++)
	{
		for(m=0;m<nodes_x;m++)
		{
			
			array_index=m+nodes_x*n;

			
			for(i=0;i<3;i++)
			{
				for(j=0;j<3;j++)
				{
					//initialising the sigma_elastic tensor
					sigma_elastic[i][j]=0.0;
					//creating the het strain tensor
					if(i==0 && j==0)
						het[i][j]=creal(het11[array_index]);
					else if(i==1 && j==1)
						het[i][j]=creal(het22[array_index]);
					else if((i==0 && j==1) || (i==1 && j==0))
						het[i][j]=creal(het12[array_index]);
					else
						het[i][j]=0.0;
					//creating the elastic strain matrix
					epsilon_elastic[i][j]=hom_strain[i][j]+het[i][j]-(theta_c[array_index]*epsilon_c[i][j]+theta_eta[array_index]*epsilon_eta[i][j]);
					//epsilon_elastic[i][j]=hom_strain[i][j]+het[i][j]-(theta_c[array_index]*epsilon_c[i][j]);
				}
			}
			if((fabs(epsilon_elastic[0][2])>1e-8) || (fabs(epsilon_elastic[1][2]>1e-8)) || (fabs(epsilon_elastic[2][2]>1e-8)))
			{
				printf("The elastic strains 13 ,23 and 33 are not zero\n");
				exit(1);
			}
			//computing the elastic moduli parameters
			lambda=lambda_mat+del_lambda*creal(fi[array_index]);
			mu=mu_mat+del_mu*creal(fi[array_index]);
			mu_prime=mu_prime_mat+del_mu_prime*creal(fi[array_index]);
			//computing the elastic stresses
			str(sigma_elastic,epsilon_elastic,kron_delta,lambda,mu,mu_prime);
			//some checks performed on the solution
	//		if(fabs(sigma_elastic[0][1])>1e-8)
	//		{
	//			printf("The elastic stresses 12 are not zero\n");
	//			exit(1);
	//		}
			if((fabs(sigma_elastic[0][2])>1e-8) || (fabs(sigma_elastic[1][2]>1e-8)) )
			{
				printf("The elastic stresses 13,23 and 33 are not zero\n");
				exit(1);
			}
			//transferring the data to the concerned arrays
			el_str_11[array_index]=sigma_elastic[0][0];
			el_str_22[array_index]=sigma_elastic[1][1];
			el_str_12[array_index]=sigma_elastic[0][1];
			//writing data to files
			//if(n==nodes_y/2 && m>=nodes_x/2)
			//	fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",((m-nodes_x/2)*dx)/radius,sigma_elastic[0][0]/(mu_mat*(epsilon_c_mag+epsilon_eta_mag)),sigma_elastic[1][1]/(mu_mat*(epsilon_c_mag+epsilon_eta_mag)),sigma_elastic[0][1]/(mu_mat*(epsilon_c_mag+epsilon_eta_mag)));

			//computing the elastic chemical potential which enters the Cahn-Hilliard equation
			mu_el_comp[array_index]=0.0+_Complex_I*0.0;
			mu_el_eta[array_index]=0.0+_Complex_I*0.0;
			//forming the del_sigma_el matrix
			str(del_sigma_el,epsilon_elastic,kron_delta,del_lambda,del_mu,del_mu_prime);
			for(i=0;i<3;i++)
			{
				for(j=0;j<3;j++)
				{	
					mu_el_comp[array_index]+=-(sigma_c[i][j]+fi[array_index]*del_sigma_c[i][j])*(1.0/(av_comp_beta-av_comp_alpha))*epsilon_elastic[i][j];
					mu_el_eta[array_index]+=-(sigma_eta[i][j]+fi[array_index]*del_sigma_eta[i][j])*(1.0/av_eta)*epsilon_elastic[i][j]+0.5*(1.0/av_eta)*del_sigma_el[i][j]*epsilon_elastic[i][j];
					
				}
			}	
			//if(t_steps==3)
			//{
			//	printf("mu_el_eta=%lf\tmu_el_eta_imag=%lf\n",creal(mu_el_eta[array_index]),cimag(mu_el_eta[array_index]));
			//	getchar();
			//}
		}		 				
		
	} 
	//fclose(fp);
}
