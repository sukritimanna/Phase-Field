//this is a routine for computing the equation of mechanical equilibrium under Higher Order Approximation

void inhom_mod_elas(fftw_complex *fi,fftw_complex *fou_fi,fftw_complex *theta_c,fftw_complex *fou_theta_c,fftw_complex *theta_eta,fftw_complex *fou_theta_eta,fftw_complex *fi_times_theta_c,fftw_complex *fi_times_theta_eta,fftw_plan plf,fftw_plan plf_fi,fftw_plan plb,double *omega_in_store,double *kx,double *ky,double lambda_mat,double mu_mat,double mu_prime_mat,double del_lambda,double del_mu,double del_mu_prime,double sigma_c[3][3],double sigma_eta[3][3],double del_sigma_c[3][3],double del_sigma_eta[3][3],double kron_delta[3][3],double hom_strain[3][3],fftw_complex *het11,fftw_complex *het22,fftw_complex *het12,fftw_complex *ux,fftw_complex *uy,int nodes_x,int nodes_y,double dx,double dy,double lambda_bar,double mu_bar,double mu_prime_bar,double toler)

{

	int m,n,i,j,k,iter=1;

	long array_index;
		
	double error=0.0,old_error,mu,mu_prime,lambda,del_sigma_E[3][3],average_eigen_stress[3][3],average_periodic_stress[3][3];	
	
	double omega_in[3][3],k_vector[3],real_uz,imag_uz;

	fftw_complex u[3],*product11,*product12,*product22,product[3][3],del_sigma_product[3][3],*old_ux,*old_uy;

	
	if((old_ux=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the old ux array\n");
		exit(1);
	}
	if((old_uy=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the old uy array\n");
		exit(1);
	}
	if((product11=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the product11 array\n");
		exit(1);
	}
	if((product22=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the product22 array\n");
		exit(1);
	}
	if((product12=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the product12 array\n");
		exit(1);
	}
	
		
	//computing the average eigen_stress (this is constant for all the iterations; see 'average_eigen_stress.c")
	avg_sigma_not(average_eigen_stress,sigma_c,sigma_eta,del_sigma_c,del_sigma_eta,theta_c,theta_eta,fi_times_theta_c,fi_times_theta_eta,nodes_x,nodes_y,dx,dy);
	//printf("The average eigen stress is\n");
	//prt(average_eigen_stress);

	//taking the fourier transform of the other interpolant
	fftw_execute_dft(plf,fi_times_theta_c,fi_times_theta_c);
	fftw_execute_dft(plf,fi_times_theta_eta,fi_times_theta_eta);
	//taking the fourier transform of the fi array
	fftw_execute(plf_fi);

	do	
	{
		//calculating the average periodic stresses (see "average_periodic_stress.c")
		avg_sigma_star(average_periodic_stress,kron_delta,het11,het12,het22,fi,lambda_mat,mu_mat,mu_prime_mat,del_lambda,del_mu,del_mu_prime,nodes_x,nodes_y,dx,dy,iter);  
		//printf("The average periodic stress at iteration %d\n",iter);
		//prt(average_periodic_stress);
	//	printf("Here\n");
	//	getchar();

		//calculating the homogeneous strain tensor (this assumes that the applied stress on the system is zero; see "homogeneous.c")

		hom_str(hom_strain,average_periodic_stress,average_eigen_stress,kron_delta,lambda_bar,mu_bar,mu_prime_bar,iter);
		//printf("The homgeneous strain at iteration=%d\n",iter); 
		//prt(hom_strain);		
		//	exit(1);
		
				
		
		//forming the stress tensor which is obtained from operating the delta modulus matrix on homogeneous strain (see "stress.c") 

		

		str(del_sigma_E, hom_strain, kron_delta, del_lambda, del_mu, del_mu_prime);
		
		//printf("The stress tensor formed by the delta modulus parameters and the homogeneous strain is\n"); 
		//prt(del_sigma_E);

		//copying the new strains into the old strains
		//and creating the product arrays
		for(n=0;n<nodes_y;n++)
		{
			for(m=0;m<nodes_x;m++)
			{
				array_index=m+nodes_x*n;

				
				old_ux[array_index]=ux[array_index];
				//printf("imag_ux=%lf\n",cimag(ux[array_index]));
				//getchar();
				old_uy[array_index]=uy[array_index];
				

				product11[array_index]=fi[array_index]*het11[array_index];
				
				product22[array_index]=fi[array_index]*het22[array_index];				
				product12[array_index]=fi[array_index]*het12[array_index];				
				//printf("product11=%lf\tproduct22=%lf\tproduct=%lf\n",cimag(product11[array_index]),cimag(product22[array_index]),cimag(product12[array_index]));
				//getchar();
				
			}
		}
		//taking the fourier transform of the product arrays	
		fftw_execute_dft(plf,product11,product11);
		fftw_execute_dft(plf,product22,product22);
		fftw_execute_dft(plf,product12,product12);
		//solving for the periodic displacements
		error=0.0;
		for(n=0;n<nodes_y;n++)
		{
			for(m=0;m<nodes_x;m++)
			{
				array_index=m+nodes_x*n;
				//creating the product matrix for this node
				for(i=0;i<3;i++)
				{
					for(j=0;j<3;j++)
					{
						if(i==0 && j==0)
							product[i][j]=product11[array_index];
						else if(i==1 && j==1)
							product[i][j]=product22[array_index];
						else if((i==0 && j==1) || (i==1 && j==0))
							product[i][j]=product12[array_index];
						else
							product[i][j]=0.0+_Complex_I*0.0;
					}
				}
				//forming the stress matrices corresponding to the product arrays (see "stress_product.c")
			//	printf("The product matrix\n");
				stress_prod(del_sigma_product,product,kron_delta,del_lambda,del_mu,del_mu_prime);
			//	getchar();
				//extracting the omega inverse matrix
				for(i=0;i<3;i++)
				{
					for(j=0;j<3;j++)
					{
						omega_in[i][j]=omega_in_store[j+(3*i)+(3*3*m)+(3*3*nodes_x*n)];
					}
				}
				//initializing the k_vector
				k_vector[0]=kx[m];
				k_vector[1]=ky[n];
				k_vector[2]=0.0;//this line could have been placed outside	
				//initializing the u_vector
				u[0]=0.0+_Complex_I*0.0;	
				u[1]=0.0+_Complex_I*0.0;	
				u[2]=0.0+_Complex_I*0.0;
				//printf("real_fou_theta_eta=%lfimag_fou_theta=%lf\n",creal(fou_theta_eta[array_index]),cimag(fou_theta_eta[array_index]));
				//printf("real_fou_theta_c=%lfimag_fou_theta_c=%lf\n",creal(fou_theta_c[array_index]),cimag(fou_theta_c[array_index]));
				//printf("real_fou_fi=%lfimag_fou_fi=%lf\n",creal(fou_fi[array_index]),cimag(fou_fi[array_index]));
				//printf("real_fi_times_theta_c=%lfimag_fi_times_theta_c=%lf\n",creal(fi_times_theta_c[array_index]),cimag(fi_times_theta_c[array_index]));
				//printf("real_fi_times_theta_eta=%lfimag_fi_times_theta_eta=%lf\n",creal(fi_times_theta_eta[array_index]),cimag(fi_times_theta_eta[array_index]));
		//		printf("real fi_times_theta_eta=%lf\timag_fi_times_theta_eta=%lf\n",creal(fi_times_theta_eta[array_index]),cimag(fi_times_theta_eta[array_index]));
		//		getchar();
				//starting the computations
				for(k=0;k<3;k++)
				{
					for(i=0;i<3;i++)
					{
						for(j=0;j<3;j++)
						{
							//we have removed the negative sign of the formula as this has already been taken care of in the omega routine
							//printf("real_del_sigma_prod=%lf\timag_del_sigma_prod=%lf\n",creal(del_sigma_product[i][j]),cimag(del_sigma_product[i][j]));
							u[k]+=_Complex_I*k_vector[j]*omega_in[k][i]*(sigma_c[i][j]*fou_theta_c[array_index]+sigma_eta[i][j]*fou_theta_eta[array_index]-del_sigma_E[i][j]*fou_fi[array_index]+del_sigma_c[i][j]*fi_times_theta_c[array_index]+del_sigma_eta[i][j]*fi_times_theta_eta[array_index]-del_sigma_product[i][j]);
							//u[k]+=_Complex_I*k_vector[j]*omega_in[k][i]*(sigma_c[i][j]*fou_theta_c[array_index]-del_sigma_E[i][j]*fou_fi[array_index]+del_sigma_c[i][j]*fi_times_theta_c[array_index]-del_sigma_product[i][j]);
						}
					}
				}
				
				//transferring the solutions to the disp array
				ux[array_index]=u[0];
				//printf("imag_u[0]=%lf\n",cimag(u[0]));
				uy[array_index]=u[1];
				
				real_uz=creal(u[2]);
				imag_uz=cimag(u[2]);
				if(fabs(real_uz)!=0.0 || fabs(imag_uz)!=0.0)
				{
					printf("The u_z displacements are not zero in the fourier space for the zeroth order approximation\n");
					exit(1);
				}
				
				error+=(creal(ux[array_index])-creal(old_ux[array_index]))*(creal(ux[array_index])-creal(old_ux[array_index]));
				//printf("real_ux=%lf\treal_old_ux=%lf\n",creal(ux[array_index]),creal(old_ux[array_index]));
				error+=(cimag(ux[array_index])-cimag(old_ux[array_index]))*(cimag(ux[array_index])-cimag(old_ux[array_index]));
				//printf("imag_ux=%lf\timag_old_ux=%lf\n",cimag(ux[array_index]),cimag(old_ux[array_index]));
				error+=(creal(uy[array_index])-creal(old_uy[array_index]))*(creal(uy[array_index])-creal(old_uy[array_index]));
				//printf("real_ux=%lf\treal_old_ux=%lf\n",creal(ux[array_index]),creal(old_ux[array_index]));
				error+=(cimag(uy[array_index])-cimag(old_uy[array_index]))*(cimag(uy[array_index])-cimag(old_uy[array_index]));
				//printf("imag_ux=%lf\timag_old_ux=%lf\n",cimag(ux[array_index]),cimag(old_ux[array_index]));
				//printf("error=%lf\n",error);
				//getchar();
				//printf("array_index=%ld\n",array_index);
			}
		}
		
		//forming the heterogeneous strain arrays
				
		het_strains(ux,uy,het11,het22,het12,plb,kx,ky,nodes_x,nodes_y,iter);	
		error=sqrt(error*dx*dy);
		//printf("error=%.12lf\n",error);
		iter++;
		if(iter>2)
		{
			if(error>old_error)
			{
				printf("solution diverges at iter=%d\n",iter);
				exit(1);
			}
		}
		old_error=error;
					
	}
	while(error>toler);
	//freeing up the arrays
	fftw_free(old_ux);
	fftw_free(old_uy);
	
	fftw_free(product11);
	fftw_free(product22);	
	fftw_free(product12);
}
