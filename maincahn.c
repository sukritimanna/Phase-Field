//this contains the 'main' of the code which computes the elastic stress fields for a circular precipitate placed at the centre of a square domain of the matrix
//we measure strains with respect to the lattice parameter of the matrix 
//again under situations of inhomogeneous modulus we measure the modulus with respect to the matrix
//we are going to consider materials with cubic anisotropy from which the isotropic situation can be derived
//the main objective of this code is the testing of the new shortened elasticity formulation
//we are going to create a provision such that there is an eigen strain contributuion from both c and eta
//we are going to use the equilibrium composition and eta of the different phases from the binary phase diagram we have
//the elastic modulus is a function of c alone
//we are 
//the constraint is that the applied stress is zero
// we are going to employ a plane strain approximation
//for the time being we will assume that the eigen strains are dilatational 

//including the header files

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<fftw3.h>
#include <time.h>
#include <float.h>
//#define p 5000 /* total no pf generated random location*/			
//#define CIRC 200  /*no of particles*/

#include"inputparamcahn.c"
#include"compini_single_particle.c"
#include"compini_single_particle_4Grain.c"
#include"compini_multi_particle.c"
#include"k_vector_gen.c"
#include"matrix_printer.c"
#include"epsilon_matrices.c"
#include"elas_const.c"
#include"matrix_inversion.c"
#include"storage.c"
#include"omega.c"
#include"heterogeneous.c"
#include"hom_mod_elas.c"
#include"average.c"
#include"compliance.c"
#include"stress.c"
#include"average_eigen_stress.c"
#include"average_periodic_stress.c"
#include"homogeneous.c"
#include"stress_product.c"
#include"inhom_mod_elas.c"
#include"elastic_chem_pot.c"
//#include"radius_cal.c"


//****************************************************************************//
main()
{
system("rm -f ../output/data/comp*");//To delete all the old files
system("rm -f ../output/data/eta*");//To delete all the old files
	char fname[300];
	char NAME[300];
	FILE *fp,*fp_elas,*fp_1d,*fp_rad;
	FILE *output4;//Binary-File Saving for 'e'
	double sys_size_x,sys_size_y,dx,dy,nu,A_z,delta,C_11_mat,C_12_mat,C_44_mat;

	double C_11_ppt,C_12_ppt,C_44_ppt,del_C_11,del_C_12,del_C_44,avg_C_11,avg_C_12,avg_C_44;

	double lambda_mat,mu_mat,mu_prime_mat,lambda_ppt,mu_ppt,mu_prime_ppt,del_lambda,del_mu,del_mu_prime;
			
	double S_11,S_12,S_44,lambda_bar,mu_bar,mu_prime_bar;

	double *kx,*ky,*sq_mag_k_vectors,epsilon_c_mag,epsilon_eta_mag,toler;

	double *omega_in_store,radius,epsilon_c[3][3],epsilon_eta[3][3],kron_delta[3][3];

	double sigma_c[3][3],sigma_eta[3][3],del_sigma_c[3][3],del_sigma_eta[3][3];
	
	double hom_strain[3][3];
	
	double M,L,K_eta,K_comp,dt;
		
	//int total_timesteps,output_timesteps_interval,;
	
	int total_time,output_time_interval,t_steps,num_eta,num,INDEX;
	
	fftw_complex *comp,*eta,*fi,*fou_fi,*theta_c,*fou_theta_c,*theta_eta,*fou_theta_eta,*fi_times_theta_c,*fi_times_theta_eta; //though there is no explicit need for the comp and e ta arrays we declare them just as a debugging aid
		
	
	double av_comp_alpha,av_comp_beta,av_eta,temp,delta_f;
		
	double comp_real,eta_real; 

	double *el_str_11,*el_str_22,*el_str_12,*c,*sum;

	int nodes_x,nodes_y,half_mark_x,half_mark_y,i,j,m,n;

	fftw_plan plf,plf_fi,plb;
	
	double realc,realeta,realc1;
	double realetam,realetan;
	double eta_pow1,eta_pow2,eta_pow3,sum2;
	double etaietaj_pow2,etapow2_sum,etapow4_sum,etasum_gb;
	double A,B,Z;
	long array_index;
        
      //  double rad,x1;
	
        fftw_complex *het11,*het22,*het12,*ux,*uy,*mu_el_eta,*mu_el_comp,*g_func,*v_func;
	//computing the values of omega_1 and omega_2
	//double omega_1,omega_2;
	//omega_1=N*z_1*e_1*ev2Joule;
	//omega_2=N*z_2*e_2*ev2Joule;
	
	//printf("omega_1=%lf\tomega_1=%lf\n",omega_1,omega_2);	

	//reading the system parameters (see "inputparamcahn.c")
	inputparam(&sys_size_x,&sys_size_y,&dx,&dy,&M,&L,&K_comp,&K_eta,&epsilon_c_mag,&epsilon_eta_mag,&C_44_mat,&nu,&A_z,&delta,&radius,&av_comp_alpha,&av_comp_beta,&av_eta,&toler,&dt,&total_time,&output_time_interval,&num_eta,&A,&B,&Z);


	//calculating the no. of nodes in the system
	nodes_x=sys_size_x/dx;
	nodes_y=sys_size_y/dy;
	printf("The no. of nodes in the system's x-direction will be=%d\n",nodes_x);
	printf("The no. of nodes in the system's y-direction will be=%d\n",nodes_y);
	//calculating the time stepping information
	//total_timesteps=(int)(total_time/dt);
	//output_timesteps_interval=(int)(output_time_interval/dt);
	//printf("The total no. of timesteps the simulation will run=%d\n",total_timesteps);
	//printf("The total no. of timesteps at which the output will be recorded=%d\n",output_timesteps_interval);


	//creating different arrays
	if((comp=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the composition array\n");
		exit(1);
	}
	if((eta=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the eta array\n");
		exit(1);
	}
	if((g_func=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the g_func array\n");
		exit(1);
	}
	if((v_func=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the v_func array\n");
		exit(1);
	}
	if((mu_el_comp=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the mu_el_comp array\n");
		exit(1);
	}
	if((mu_el_eta=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the mu_el_eta array\n");
		exit(1);
	}
	if((fi=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the fi array\n");
		exit(1);
	}						
	if((fou_fi=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the fou_fi array\n");
		exit(1);
	}						
	if((theta_c=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the theta_c array\n");
		exit(1);
	}
	if((fou_theta_c=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the fou_theta_c array\n");
		exit(1);
	}	
	if((theta_eta=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the theta_eta array\n");
		exit(1);
	}
	if((fou_theta_eta=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the fou_theta_eta array\n");
		exit(1);
	}
	if((fi_times_theta_eta=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the fi_times_theta_eta array\n");
		exit(1);
	}
	if((fi_times_theta_c=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the fi_times_theta_c array\n");
		exit(1);
	}
	if((het11=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the het11 array\n");
		exit(1);
	}
	if((het22=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the het22 array\n");
		exit(1);
	}
	if((het12=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the het12 array\n");
		exit(1);
	}
	if((ux=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the ux array\n");
		exit(1);
	}
	if((uy=(fftw_complex *)fftw_malloc(nodes_x*nodes_y*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the uy array\n");
		exit(1);
	}
	
	if((kx=(double *)malloc(nodes_x*sizeof(double)))==NULL)
	{
		printf("Space creation has failed for the kx vectors array\n");
		exit(1);
	}
	if((ky=(double *)malloc(nodes_y*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the ky vectors array\n");
		exit(1); 
	}
	if((sq_mag_k_vectors=(double *)malloc(nodes_x*nodes_y*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the squared magnitude of the k vectors array\n");
		exit(1);
	}
	if((el_str_11=(double *)malloc(nodes_x*nodes_y*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the el_str_11 array\n");
		exit(1);
	}
	if((el_str_22=(double *)malloc(nodes_x*nodes_y*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the el_str_22 array\n");
		exit(1);
	}
	if((el_str_12=(double *)malloc(nodes_x*nodes_y*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the el_str_12 array\n");
		exit(1);
	}
		
	if((omega_in_store=(double *)malloc(3*3*nodes_x*nodes_y*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the stoarge array\n");
		exit(1);
	}
	if((c=(double *)malloc(nodes_x*nodes_y*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the c array\n");
		exit(1);
	}
	if((sum=(double *)malloc(nodes_x*nodes_y*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the sum array\n");
		exit(1);
	}
	fftw_complex *eta1[num_eta];
	fftw_complex *g_eta[num_eta];
	for(num=0;num<num_eta;num++)
	{
	eta1[num] = fftw_malloc(nodes_x*nodes_y*num_eta*sizeof(fftw_complex));
	g_eta[num] = fftw_malloc(nodes_x*nodes_y*num_eta*sizeof(fftw_complex));
	}
	
      /* if((comp_1d=(fftw_complex *)fftw_malloc(nodes_x*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the composition_1_D array\n");
		exit(1);
	}
        if((eta_1d=(fftw_complex *)fftw_malloc(nodes_x*sizeof(fftw_complex)))==NULL)
	{
		printf("Space creation failed for the eta_int array\n");
		exit(1);
	}
        if((fp_rad=fopen("../DATA/radius.dat","w"))==NULL)
	{
		printf("file for radius can't be opened\n");
		exit(1);
	}*/

        fftw_init_threads();
        fftw_plan_with_nthreads(8);
	//creating the plans for the fourier transform
	plf=fftw_plan_dft_2d(nodes_y,nodes_x,comp,comp,FFTW_FORWARD,FFTW_EXHAUSTIVE);
	plb=fftw_plan_dft_2d(nodes_y,nodes_x,comp,comp,FFTW_BACKWARD,FFTW_EXHAUSTIVE);
	//creating an out of place forward for the fi array
	plf_fi=fftw_plan_dft_2d(nodes_y,nodes_x,fi,fou_fi,FFTW_FORWARD,FFTW_EXHAUSTIVE);


	//creating the initial composition profile (see "compinicahn.c")
	//compini(nodes_x,nodes_y,dx,dy,epsilon_c_mag,epsilon_eta_mag,C_44_mat,nu,A_z,delta,radius,av_comp_alpha,av_comp_beta,av_eta,comp,eta,fi,theta_c,theta_eta,fi_times_theta_c,fi_times_theta_eta,g_func,v_func,num_eta,eta1);
	compini_multi(nodes_x,nodes_y,dx,dy,epsilon_c_mag,epsilon_eta_mag,C_44_mat,nu,A_z,delta,radius,av_comp_alpha,av_comp_beta,av_eta,comp,eta,fi,theta_c,theta_eta,fi_times_theta_c,fi_times_theta_eta,g_func,v_func,num_eta,eta1);
	//compini_4Grain(nodes_x,nodes_y,dx,dy,epsilon_c_mag,epsilon_eta_mag,C_44_mat,nu,A_z,delta,radius,av_comp_alpha,av_comp_beta,av_eta,comp,eta,fi,theta_c,theta_eta,fi_times_theta_c,fi_times_theta_eta,g_func,v_func,num_eta,eta1);
	//we have initialised the interpolant arrays in the 'compini' function itself so that we can directly solve the homogeneous modulus elasticity problem before entering the tome loop
	

	//creating the k-vector arrays (see "k_vector_gen.c")
	k_vec(kx, ky, sq_mag_k_vectors, nodes_x,nodes_y, sys_size_x, sys_size_y); 


	//creating the epsilon_c, epsilon_eta and kron_delta matrices (see "epsilon_matrices.c")
	epsilon_mat(epsilon_c,epsilon_eta, kron_delta, epsilon_c_mag,epsilon_eta_mag);

	//printing the different elasticity related matrices (see "matrix_printer.c")
	printf("The epsilon_c matrix is \n");
	prt(epsilon_c);

	printf("The epsilon_eta matrix is \n");
	prt(epsilon_eta);

	printf("The kron_delta matrix is \n");
	prt(kron_delta);
	
	//creating the elastic constants (see "elas_const.c")
	el_const( C_44_mat,delta, nu, A_z,&C_12_mat, &C_11_mat,&C_12_ppt, &C_11_ppt,&C_44_ppt,&lambda_mat,&mu_mat,&mu_prime_mat,&lambda_ppt,&mu_ppt,&mu_prime_ppt,&del_lambda,&del_mu,&del_mu_prime,&del_C_12, &del_C_11,&del_C_44); 
	
		
	//printing the elastic constants
	//printf("The elastic constants for the matrix is:C_12=%lf\tC_11=%lf\tC_44=%lf\nlambda=%lf\tmu=%lf\tmu_prime=%lf\n",C_12_mat,C_11_mat,C_44_mat,lambda_mat,mu_mat,mu_prime_mat);
	
	//printf("The elastic constants for the precipitate phase is:C_12=%lf\tC_11=%lf\tC_44=%lf\nlambda=%lf\tmu=%lf\tmu_prime=%lf\n",C_12_ppt,C_11_ppt,C_44_ppt,lambda_ppt,mu_ppt,mu_prime_ppt);

	//printf("The difference in the modulus related parameters\n");
	
	//printf("del_lambda=%lf\tdel_mu=%lf\tdel_mu_prime=%lf\n",del_lambda,del_mu,del_mu_prime); 
	
	//printf("del_C_12=%lf\tdel_C_11=%lf\tdel_C_44=%lf\n",del_C_12,del_C_12,del_C_44); 

          

	//creating the omega and omega_in matrices and storing them up (see "omega.c")
 	omega_cr(omega_in_store,nodes_x,nodes_y,kx,ky,sq_mag_k_vectors,lambda_mat,mu_mat,mu_prime_mat,kron_delta);	
        //printf("code_correct_2");
	//computing the different sigma tensors (see "stress.c")

	printf("The sigma_c tensor is\n");
	str(sigma_c,epsilon_c,kron_delta,lambda_mat,mu_mat,mu_prime_mat);
	prt(sigma_c);

	printf("The sigma_eta tensor is\n");
	str(sigma_eta,epsilon_eta,kron_delta,lambda_mat,mu_mat,mu_prime_mat);
	prt(sigma_eta);

	printf("The del_sigma_c tensor is\n");
	str(del_sigma_c,epsilon_c,kron_delta,del_lambda,del_mu,del_mu_prime);
	prt(del_sigma_c);

	printf("The del_sigma_eta tensor is\n");
	str(del_sigma_eta,epsilon_eta,kron_delta,del_lambda,del_mu,del_mu_prime);
	prt(del_sigma_eta);
	
	//computing the homogeneous elasticity problem (see "hom_mod_elas.c")
	hom_elas(theta_c,fou_theta_c,theta_eta,fou_theta_eta,plf_fi,plb,sigma_c,sigma_eta,kx,ky,omega_in_store,het11,het22,het12,ux,uy,nodes_x,nodes_y);
       // printf("code_correct_1");
       // getchar();
	
//*******************************************************************************************************//
//the above partition demarcates the operations which are set once and do not change throughout the run from those which have to be updated at every timestep, in the context of a phase field code   	
	//starting the time loop
	for(t_steps=0;t_steps<=total_time;t_steps++)
	{
                 printf("%d\n",t_steps);
				 if(t_steps%1==0)
				{
					
					printf ("Time=%ld\n",t_steps);
					printf ("%lf\n",creal(comp[nodes_x/2]));
					output4=fopen("../output/data/comp_status.dat","a");
					fprintf(output4,"\nTime=%ld\tComp=%lf\n",t_steps,creal(comp[nodes_x/2]));
					fclose(output4);
				}
                // getchar();
		//THE ELASTICITY PART
		//**********************************************************************************
		//solving for the heterogeneous strains
		//we have to compute the average C_11,C_12 and C_44 over the entire domain (see "average.c") 
		
		avg(&avg_C_11,&avg_C_12,&avg_C_44,nodes_x,nodes_y,dx,dy,C_44_mat,del_C_44,nu,A_z,fi);
			
		//creating the compliance tensor components S_11,S_12 and S_44 from the avg_C_11,avg_C_12 and avg_C_44 (see "compliance.c")

		compl(&S_11,&S_12,&S_44,&lambda_bar,&mu_bar,&mu_prime_bar,avg_C_11,avg_C_12,avg_C_44);

		//we are know going to solve the equation of mechanical equilibrium using Higher order approximations(see "inhom_mod_elas.c" )

		inhom_mod_elas(fi,fou_fi,theta_c,fou_theta_c,theta_eta,fou_theta_eta,fi_times_theta_c,fi_times_theta_eta,plf,plf_fi,plb,omega_in_store,kx,ky, lambda_mat,mu_mat,mu_prime_mat,del_lambda,del_mu,del_mu_prime,sigma_c,sigma_eta,del_sigma_c,del_sigma_eta, kron_delta, hom_strain,het11,het22,het12,ux,uy,nodes_x,nodes_y,dx,dy,lambda_bar,mu_bar,mu_prime_bar,toler);

		//the calculation of elastic chemical potentials

		 elas_chem_pot(het11,het22,het12,fi,theta_c,theta_eta,hom_strain,kron_delta,epsilon_c, epsilon_eta,epsilon_c_mag,epsilon_eta_mag,lambda_mat,mu_mat,mu_prime_mat,del_lambda,del_mu, del_mu_prime,  nodes_x, nodes_y, A_z, delta,radius,dx,mu_el_comp,mu_el_eta,el_str_11,el_str_22,el_str_12,sigma_c,sigma_eta,del_sigma_c,del_sigma_eta,av_comp_alpha,av_comp_beta,av_eta,t_steps);
		// we have some redundant arguments like"radius","dx","epsilon_c_mag",epsilon_eta_mag" 
		//in the above function call (as we have created the above function by modifying "elastic_stress.c" from the one shot calculations
		
		//************************************************************************	
		//
		//g and v func calculation
		for(j=0;j<nodes_y;j++)
		{
			for(i=0;i<nodes_x;i++)
			{
				array_index=i+nodes_x*j;

			
				comp_real = creal(comp[array_index]);
				etaietaj_pow2=0.0;
				etapow2_sum=0.0;
				etapow4_sum=0.0;
				etasum_gb=0.0;
				
            	for(m=0;m<num_eta;m++)
				{
                	for(n=m+1;n<num_eta;n++)
                    		{
				     	realetam=creal(eta1[m][array_index]);
					 	realetan=creal(eta1[n][array_index]);
                     	etaietaj_pow2=etaietaj_pow2+realetam*realetam*realetan*realetan;
                    		}
							
					realeta=creal(eta1[m][array_index]);
                	etapow2_sum = etapow2_sum+realeta*realeta;
                	etapow4_sum = etapow4_sum+realeta*realeta*realeta*realeta;
                	etasum_gb = etasum_gb + realeta*(1.0-realeta)*realeta*(1.0-realeta);
                }
				
				/*calculating g*/
		 		g_func[array_index]=A*2*comp_real*(1-comp_real)*(1-2*comp_real)-B*2*(1.0-comp_real)*(0.25*(etapow4_sum)-0.5*(etapow2_sum)+2*etaietaj_pow2+0.25)+2*Z*comp_real*etapow2_sum+_Complex_I*0.0;

			/*calculating g_eta*/
				for(num=0;num<num_eta;num++)
       				{
           			eta_pow1=creal(eta1[num][array_index]);
           			eta_pow2=eta_pow1*eta_pow1;
           			eta_pow3=eta_pow2*eta_pow1;
         			g_eta[num][array_index]=B*(1.0-comp_real)*(1.0-comp_real)*(eta_pow3-eta_pow1+4.0*eta_pow1*(etapow2_sum-eta_pow2))+2*Z*eta_pow1*comp_real*comp_real+_Complex_I*0.0;
					}
			}
			
		}

		//opening the files for writing the composition and eta fields in one and the stress fields in  another
		if(t_steps%output_time_interval==0)
		{
			sprintf(fname,"../DATA_11/comp_and_eta_ep_c_%.4lfep_eta_%.4lf_A_z_%.1lf_delta_%.1lf_timestep_%d.dat",epsilon_c_mag,epsilon_eta_mag,A_z,delta,t_steps);
			if((fp=fopen(fname,"w"))==NULL)
			{
				printf("comp and eta profile can't be written\n");
				exit(1);
			}
			printf ("Time=%ld Writing to output/data/ \n",t_steps);
					output4=fopen("../output/data/comp_status.dat","a");
					fprintf(output4,"\nTime=%ld\t Writing Files to output/data/\n",t_steps);
					fclose(output4);
/*-------------------------Writing density(start)----------------------*/
		FILE *output3;//Binary-File Saving for 'e'
		for(j=0;j<nodes_y;j++)
		{
			for(i=0;i<nodes_x;i++)
			{
				array_index=i+nodes_x*j;
				sum[array_index]=0.0;
				for (INDEX=0;INDEX<=num_eta-1;INDEX++)
				{
					sum[array_index]=sum[array_index]+(INDEX+25)*pow(creal(eta1[INDEX][array_index]),2);
				}
				sum[array_index]=(sum[array_index]+(1.0-creal(comp[array_index])))/(1.0*num_eta);
				
			}
		}
		sprintf(NAME,"../output/data/comp%d.dat",(int) (t_steps));
		output3=fopen(NAME,"w");
		for(j=0;j<nodes_y;j++)
		{
			for(i=0;i<nodes_x;i++)
			{
				array_index=i+nodes_x*j;

				(c[array_index]) = creal(comp[array_index]);

				}
		}
		fwrite(&c[0],sizeof(double),(size_t)nodes_x*nodes_y,output3);
		fclose(output3);
		FILE *output2;//Binary-File Saving for 'e'
		sprintf(NAME,"../output/data/eta%d.dat",(int) (t_steps));
		output2=fopen(NAME,"w");
		fwrite(&sum[0],sizeof(double),(size_t)nodes_x*nodes_y,output2);
		fclose(output2);
		for(j=0;j<nodes_y;j++)
		{
			for(i=0;i<nodes_x;i++)
			{
				array_index=i+nodes_x*j;

				(c[array_index]) = creal(eta1[0][array_index]);

				}
		}
		FILE *output4;//Binary-File Saving for 'e'
		sprintf(NAME,"../output/data/etaC%d.dat",(int) (t_steps));
		output4=fopen(NAME,"w");
		fwrite(&c[0],sizeof(double),(size_t)nodes_x*nodes_y,output2);
		fclose(output4);
		//printf ("a\n");
		}
		//printf ("b\n");
		//THE CAHN_HILLIARD AND ALLEN_CAHN PART
		//************************************************************************************
		//moving the comp and eta arrays into fourier space
		fftw_execute_dft(plf,comp,comp);
		fftw_execute_dft(plf,eta,eta);
		//moving the g_func, v_func, mu_el_comp and mu_el_eta into the fourier space
		fftw_execute_dft(plf,mu_el_comp,mu_el_comp);		 
		fftw_execute_dft(plf,mu_el_eta,mu_el_eta);
		fftw_execute_dft(plf,g_func,g_func);
		fftw_execute_dft(plf,v_func,v_func);
		for(num=0;num<num_eta;num++)
		{
		fftw_execute_dft(plf,eta1[num],eta1[num]);
		fftw_execute_dft(plf,g_eta[num],g_eta[num]);
		}
		//printf ("c\n");
		//starting the space loops and solving the Cahn-Hilliard and the Allen-Cahn     equations		
		for(j=0;j<nodes_y;j++)
		{
			for(i=0;i<nodes_x;i++)
			{
				array_index=i+nodes_x*j;


				//if(t_steps==5)
				//{
				//	printf("fi=%lf+I*(%lf)",creal(mu_el_comp[array_index]),cimag(mu_el_comp[array_index]));
		//			printf("fi=%lf+I*(%lf)",creal(mu_el_eta[array_index]),cimag(mu_el_eta[array_index]));
		//			getchar();

		//		}
				//printf("%ld\n",array_index);
				//getchar();
				if (creal(pow((eta1[0][array_index]),2.0))>0.01 && creal(pow((eta1[0][array_index]),2.0)) <0.99)
				{
					mu_el_comp[array_index]=0.0+_Complex_I*0.0;
				}
				comp[array_index]=(comp[array_index]-(sq_mag_k_vectors[array_index]*M*dt*(g_func[array_index]+mu_el_comp[array_index])))/(1.0+(2.0*M*K_comp*dt*sq_mag_k_vectors[array_index]*sq_mag_k_vectors[array_index]));			
				if (creal(comp[array_index])<0.01 && creal(comp[array_index]) >0.99)
				{
					mu_el_eta[array_index]=0.0+_Complex_I*0.0;
				}
				for(num=0;num<num_eta;num++)
				{
				eta1[num][array_index]=(eta1[num][array_index]-(L*dt*(g_eta[num][array_index]+mu_el_eta[array_index])))/(1.0+(2.0*L*K_eta*sq_mag_k_vectors[array_index]*dt));
				}
				//if(t_steps==1)
				//{
				//	printf("fi=%lf+I*(%lf)",creal(comp[array_index]),cimag(comp[array_index]));
		//			printf("fi=%lf+I*(%lf)",creal(eta[array_index]),cimag(eta[array_index]));
		//			getchar();

		//		}
			}
		}
		//taking the inverse fourier transform of the comp and eta arrays
		fftw_execute_dft(plb,comp,comp);
		fftw_execute_dft(plb,eta,eta);
		for(num=0;num<num_eta;num++)
		{
			fftw_execute_dft(plb,eta1[num],eta1[num]);
		}
		//printf ("c\n");
		//normalising the comp and eta arrays and creating the different interpolants/functions for the next timestep (also writing the files for current time step)
		for(j=0;j<nodes_y;j++)
		{
			for(i=0;i<nodes_x;i++)
			{	
				array_index=i+nodes_x*j;
				//normalising the arrays
				comp[array_index]=comp[array_index]/(double)(nodes_x*nodes_y);		
				eta[array_index]=eta[array_index]/(double)(nodes_x*nodes_y);
				for(num=0;num<num_eta;num++)
				{
				eta1[num][array_index]=eta1[num][array_index]/(double)(nodes_x*nodes_y);
				}
				//writing the files
				if(t_steps%output_time_interval==0)
				{
					fprintf(fp,"%d\t%d\t%lf\t%lf\n",i,j,creal(comp[array_index]),creal(eta[array_index]));
					//fprintf(fp_elas,"%d\t%d\t%lf\t%lf\t%lf\n",i,j,el_str_11[array_index]/(mu_mat*(epsilon_c_mag+epsilon_eta_mag)),el_str_22[array_index]/(mu_mat*(epsilon_c_mag+epsilon_eta_mag)),el_str_12[array_index]/(mu_mat*(epsilon_c_mag+epsilon_eta_mag)));
				}						
                                  


				//creating the interpolants
				theta_c[array_index]=(comp[array_index]); //this is defined in a way such that it takes a value of 1 at the beta ppt. composition and a value of zero in tha matrix
				theta_eta[array_index]=pow((eta1[0][array_index]),2.0); //this is defined in a way such that it takes a value of 1 at the beta ppt. eta and a value of zero in tha matrix
				fi[array_index]=theta_eta[array_index]; //so fi is a function of eta
				//if(t_steps==1)
				//{
					//printf("fi=%lf+I*(%lf)",creal(fi[array_index]),cimag(fi[array_index]));
					//getchar();

				//}
				fi_times_theta_eta[array_index]=fi[array_index]*theta_eta[array_index];

				fi_times_theta_c[array_index]=fi[array_index]*theta_c[array_index];
				//creating the g_function and the v_function
				comp_real=creal(comp[array_index]);
				eta_real=creal(eta[array_index]);
				//g_func[array_index]= 2*comp_real*(1- 10*pow(eta_real,3.0) +15*pow(eta_real,4.0) -6*pow(eta_real,5.0)) -2*(1-comp_real)*(10*pow(eta_real,3.0) -15*pow(eta_real,4.0) +6*pow(eta_real,5.0)) + _Complex_I*0.0;
				//v_func[array_index]= -30*comp_real*comp_real*(pow(eta_real,2.0)- 2*pow(eta_real,3.0) + pow(eta_real,4.0)) +30*(1-comp_real)*(1-comp_real)*(pow(eta_real,2.0)- 2*pow(eta_real,3.0) + pow(eta_real,4.0)) + 2*eta_real*(1- eta_real)*(1-eta_real) - 2*eta_real*eta_real*(1-eta_real) + _Complex_I*0.0;
			//	if(t_steps==5)
			//	{
			//		printf("fi=%lf+I*(%lf)",creal(g_func[array_index]),cimag(g_func[array_index]));
		//			printf("fi=%lf+I*(%lf)",creal(v_func[array_index]),cimag(v_func[array_index]));
		//			getchar();

//				}
		
			}
			if(t_steps%output_time_interval==0)
			{
				fprintf(fp,"\n");
				//fprintf(fp_elas,"\n");
			}
		}
                        //writing data for 1 D profile
                              /*  if(t_steps%output_time_interval==0)
                                  {
                                    for(i=0;i<nodes_x;i++)
	                                   {
		                            comp_1d[i]=comp[i+nodes_x*(nodes_y/2)];
		                            eta_1d[i]=eta[i+nodes_x*(nodes_y/2)];
	                                    fprintf(fp_1d,"%d\t%lf\t%lf\n",i,creal(comp_1d[i]),creal(eta_1d[i]));
                                            }

                                  }
                              if (t_steps%output_time_interval==0)
                                 {
                                        for  (i=0;i<nodes_x/2;i++)
                                   { 
              
                                      if   (creal(comp_1d[i])<0.5&&creal(comp_1d[i+1])>0.5)
                                    {
                                  //  printf("%lf\t%lf\n%d\t%d\n",comp_int[i],comp_int[i+1],i,i+1);
		                    //x1 = (i+1)- ((comp_int[i+1]-0.65)/(comp_int[i+1]-comp_int[i]));
                                      x1=(2*i+1)*0.5;
		                     // printf("%lf\n",x1);
	
                                    	break;
                                      }
                                      }
                                    rad=(nodes_x*0.5-x1)*0.5;
                                    fprintf(fp_rad,"%d\t%lf\n",t_steps,rad);
                                  }*/
		//performing fourier transforms of theta_c and theta_eta
		fftw_execute_dft(plf_fi,theta_c,fou_theta_c);
		fftw_execute_dft(plf_fi,theta_eta,fou_theta_eta);
		if(t_steps%output_time_interval==0)
		{
			fclose(fp);
			//fclose(fp_elas);
                       // fclose(fp_1d);
		}
		//getchar();
		//closing the time loop

	}			
	//fclose(fp_rad);
	//destroying the plans
	fftw_destroy_plan(plf);
	fftw_destroy_plan(plb);
        fftw_cleanup_threads();
        fftw_cleanup();
	//freeing up the arrays
	
	fftw_free(comp);	
	fftw_free(eta);
	fftw_free(g_func);	
	fftw_free(v_func);		
	fftw_free(mu_el_comp);	
	fftw_free(mu_el_eta);	
	fftw_free(fi);	
	fftw_free(fou_fi);	
	fftw_free(theta_c);
	fftw_free(fou_theta_c);	
	fftw_free(theta_eta);
	fftw_free(fou_theta_eta);	
	fftw_free(fi_times_theta_c);	
	fftw_free(fi_times_theta_eta);	
	fftw_free(het11);	
	fftw_free(het22);
	fftw_free(het12);
	for(num=0;num<num_eta;num++)
	{
	fftw_free(g_eta[num]);
	fftw_free(eta1[num]);
	}
        
       // fftw_free(comp_1d);	
	//fftw_free(eta_1d);		

	fftw_free(ux);
	fftw_free(uy);	
	
	free(kx);	
	free(ky);	
	free(sq_mag_k_vectors);	
	free(omega_in_store);
	free(el_str_11);
	free(el_str_22);
	free(el_str_12);
        
}
