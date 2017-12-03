//this is an array to create the initial composition profile 

void compini(int nodes_x,int nodes_y,double dx,double dy,double epsilon_c_mag,double epsilon_eta_mag,double C_44_mat,double nu,double A_z,double delta,double radius,double av_comp_alpha,double av_comp_beta,double av_eta,fftw_complex *comp,fftw_complex *eta,fftw_complex *fi,fftw_complex *theta_c,fftw_complex *theta_eta,fftw_complex *fi_times_theta_c,fftw_complex *fi_times_theta_eta,fftw_complex *g_func,fftw_complex *v_func,int num_eta,fftw_complex *eta1[])
{
	printf("Single Particle Code Running\n");
	srand(time(0)) ;
	int i,j,num;
	long array_index;
	double distance,comp_real,eta_real,sum2;
	FILE *fp;
	char fname[300];
	sprintf(fname,"../DATA_11/initial_profile_ep_c_%.4lfep_eta_%.4lf_A_z_%.1lf_delta_%.1lf.dat",epsilon_c_mag,epsilon_eta_mag,A_z,delta);
	if((fp=fopen(fname,"w"))==NULL)
	{
		printf("Initial profile can't be written\n");
		exit(1);
	}
	for(j=0;j<nodes_y;j++)
	{
		for(i=0;i<nodes_x;i++)
		{
			array_index=i+nodes_x*j;
			//we are creating a precipitate at the centre of the simulation domain
			for(num=0;num<num_eta;num++)
			{	
				
				
						
						sum2=(double)(rand()/(RAND_MAX*1.0));
						//eta1[num][array_index] = 0.01*(0.5-sum2);//Initiation of 'e'
						eta1[num][array_index] = 0.0+_Complex_I*0.0;//Initiation of 'e'
						
					
			}
            if (j<nodes_y/2)
			{
				eta1[0][array_index]=1.0+_Complex_I*0.0;
			}
			else
			{
				eta1[1][array_index]=1.0+_Complex_I*0.0;
			}
			distance=sqrt((i-nodes_x/2)*(i-nodes_x/2)*dx*dx + (j-nodes_y/2)*(j-nodes_y/2)*dy*dy);
                         
			if(distance<=radius)
			{
				comp[array_index]=1.0+_Complex_I*0.0;
				eta[array_index]=av_eta+_Complex_I*0.0;
				for(num=0;num<num_eta;num++)
				{	
				eta1[num][array_index]=0.0+_Complex_I*0.0;
				}
				//if (i==253&&j==nodes_y/2)
                              // {
                                 // printf("calculated_distance=%.10lf\n",distance);
                                 // printf("i=%d\tj=%d\t\n",i,j);
                                //  getchar();
                              //  }
			}
			else
			{
				comp[array_index]=0.1+_Complex_I*0.0;
				eta[array_index]=0.0+_Complex_I*0.0;
			}
			
			//theta_c[array_index]=(comp[array_index]-av_comp_alpha)/(av_comp_beta-av_comp_alpha); //this is defined in a way such that it takes a value of 1 at the beta ppt. composition and a value of zero in tha matrix
			theta_c[array_index]=comp[array_index];
			//theta_eta[array_index]=(eta[array_index])/av_eta; //this is defined in a way such that it takes a value of 1 at the beta ppt. eta and a value of zero in tha matrix
			theta_eta[array_index]=pow((eta1[0][array_index]),2.0);
			fi[array_index]=theta_eta[array_index]; //so fi is a function of eta

			fi_times_theta_eta[array_index]=fi[array_index]*theta_eta[array_index];

			fi_times_theta_c[array_index]=fi[array_index]*theta_c[array_index];
			comp_real=creal(comp[array_index]);
			eta_real=creal(eta[array_index]);
			//g_func[array_index]=2*comp_real*(1- 10*pow(eta_real,3.0) +15*pow(eta_real,4.0) -6*pow(eta_real,5.0)) -2*(1-comp_real)*(10*pow(eta_real,3.0) -15*pow(eta_real,4.0) +6*pow(eta_real,5.0)) + _Complex_I*0.0;
			//v_func[array_index]=-30*comp_real*comp_real*(pow(eta_real,2.0)- 2*pow(eta_real,3.0) + pow(eta_real,4.0)) +30*(1-comp_real)*(1-comp_real)*(pow(eta_real,2.0)- 2*pow(eta_real,3.0) + pow(eta_real,4.0)) + 2*eta_real*(1- eta_real)*(1-eta_real) - 2*eta_real*eta_real*(1-eta_real) + _Complex_I*0.0;
                        //fprintf(fp,"%d\t%d\t%lf\t%lf\t%lf\t%lf\n",i,j,creal(comp[array_index]),creal(eta[array_index]),creal(g_func[array_index]),creal(v_func[array_index]));
			fprintf(fp,"%d\t%d\t%lf\t%lf\n",i,j,creal(comp[array_index]),creal(eta[array_index]));
		}
		fprintf(fp,"\n");		
	}
	fclose(fp);
}
