//this is a rouitne for solving the equation of mechanical equilibrium under a homogeneous modulus approximation
void hom_elas(fftw_complex *theta_c,fftw_complex *fou_theta_c,fftw_complex *theta_eta,fftw_complex *fou_theta_eta,fftw_plan plf_fi, fftw_plan plb,double sigma_c[3][3],double sigma_eta[3][3],double *kx,double *ky, double *omega_in_store,fftw_complex *het11,fftw_complex *het22, fftw_complex *het12,fftw_complex *ux,fftw_complex *uy,int nodes_x,int nodes_y)
{
	long array_index;
	int m,n,i,j,k,iter=0;//an iter=0 corresponds to the homogeneous modulus approximation; we will directly set iter=1 for the full blown elasticity problem 
	double omega_in[3][3],k_vector[3],real_uz,imag_uz;
	fftw_complex u[3];
	

	//performing fourier transforms of theta_c and theta_eta
	fftw_execute_dft(plf_fi,theta_c,fou_theta_c);
	fftw_execute_dft(plf_fi,theta_eta,fou_theta_eta);
	//starting the space loops
	for(n=0;n<nodes_y;n++)
	{
		for(m=0;m<nodes_x;m++)
		{  
			array_index=m+nodes_x*n;
	//		printf("real_theta_c=%lf\timag_theta_c=%lf\n",creal(theta_c[array_index]),cimag(theta_c[array_index]));
	//		printf("real_theta_eta=%lf\timag_theta_eta=%lf\n",creal(theta_eta[array_index]),cimag(theta_eta[array_index]));
	//		printf("The omega_in tensor\n");
			//extracting the omega_inverse for this location
			for(i=0;i<3;i++)
			{
				for(j=0;j<3;j++)
				{
					omega_in[i][j]=omega_in_store[j+(3*i)+(3*3*m)+(3*3*nodes_x*n)];
	//				printf("%lf\t",omega_in[i][j]);
				}
	//			printf("\n");
			}
	//		getchar();
			//initializing the k_vector
			k_vector[0]=kx[m];
			k_vector[1]=ky[n];
			k_vector[2]=0.0;//this line could have been placed outside	
			//initializing the u_vector
			u[0]=0.0+_Complex_I*0.0;	
			u[1]=0.0+_Complex_I*0.0;	
			u[2]=0.0+_Complex_I*0.0;
			//starting the computations
			for(k=0;k<3;k++)
			{
				for(i=0;i<3;i++)
				{
					for(j=0;j<3;j++)
					{
						//we have removed the negative sign of the formula as this has already been taken care of in the omega routine
						u[k]+=_Complex_I*k_vector[j]*omega_in[k][i]*(sigma_c[i][j]*fou_theta_c[array_index]+sigma_eta[i][j]*fou_theta_eta[array_index]);
						//u[k]+=_Complex_I*k_vector[j]*omega_in[k][i]*(sigma_c[i][j]*fou_theta_c[array_index]);
					}
				}
			}
			//transferring the solutions to the disp array
			ux[array_index]=u[0];
			uy[array_index]=u[1];
			//printf("imag_ux=%lf\n",cimag(ux[array_index]));
				//getchar();
			real_uz=creal(u[2]);
			imag_uz=cimag(u[2]);
			if(fabs(real_uz)!=0.0 || fabs(imag_uz)!=0.0)
			{
				printf("The u_z displacements are not zero in the fourier space for the zeroth order approximation\n");
				exit(1);
			}
		}

	}
	//forming the heterogeneous or periodic strain arrays (see "heterogeneous.c")
	het_strains(ux,uy,het11,het22,het12,plb,kx,ky,nodes_x,nodes_y,iter);
}
