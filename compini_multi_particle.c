//this is an array to create the initial composition profile 

void compini_multi(int nodes_x,int nodes_y,double dx,double dy,double epsilon_c_mag,double epsilon_eta_mag,double C_44_mat,double nu,double A_z,double delta,double radius,double av_comp_alpha,double av_comp_beta,double av_eta,fftw_complex *comp,fftw_complex *eta,fftw_complex *fi,fftw_complex *theta_c,fftw_complex *theta_eta,fftw_complex *fi_times_theta_c,fftw_complex *fi_times_theta_eta,fftw_complex *g_func,fftw_complex *v_func,int num_eta,fftw_complex *eta1[])
{
	printf("Multi Particle Code Running\n");
	srand(time(0)) ;
	int i,j,num,a,i1,i2,IJK;
	long array_index;
	double distance,comp_real,eta_real,sum2,vf,ecc,R;
	FILE *fp;
	char fname[300];
	sprintf(fname,"../DATA_11/initial_profile_ep_c_%.4lfep_eta_%.4lf_A_z_%.1lf_delta_%.1lf.dat",epsilon_c_mag,epsilon_eta_mag,A_z,delta);
	if((fp=fopen(fname,"w"))==NULL)
	{
		printf("Initial profile can't be written\n");
		exit(1);
	}
	FILE *input;	
	if( (input = fopen("../INPAR/multi_particle.dat","r")) == NULL){
		printf("Unable to open ../INPAR/multi_particle.dat.");
		printf("Exiting\n");
		exit(0);
	}
	else{
		input=fopen("../INPAR/multi_particle.dat","r");
		printf("Reading ../INPAR/multi_particle.dat\n");
	}

	fscanf(input,"%ld %lf %lf",&a,&vf,&ecc);//Inputting Values
	fclose(input);
	printf("Reading ../INPAR/multi_particle.dat.");
double r0;
long int max_particle;
r0=1.0*a/(2);
max_particle=(int)(1.0*vf*nodes_x*nodes_y/(1.0*M_PI*r0*r0*ecc));//.................arbitrary assigning a large number (taking 1.0 as particle packing fraction i.e. not possible in reality)
printf("%ld",max_particle);
double *x,*y,*z,*rd,x1,y1,z1;
double X,Y,Z;

int INDEX;

	
x = (double *) malloc((size_t) max_particle*sizeof(double ));//Memory Allocation
y = (double  *) malloc((size_t) max_particle*sizeof(double ));
z = (double *) malloc((size_t) max_particle*sizeof(double ));
rd = (double  *) malloc((size_t)max_particle*sizeof(double ));

int nn,IN,iter=0;
double n_x_by2,n_y_by2,xxx,xxx1,xxx2,xxx3,xxx4,tmp1,tmp2,comp1,HG;

n_x_by2 = 0.5*nodes_x;	
n_y_by2 = 0.5*nodes_y;
for(j=0;j<nodes_y;j++)
	{
		for(i=0;i<nodes_x;i++)
		{
			array_index=i+nodes_x*j;
			//we are creating a precipitate at the centre of the simulation domain
			for(num=0;num<num_eta;num++)
			{	
			
						sum2=(double)(rand()/(RAND_MAX*1.0));
						eta1[num][array_index] = 0.01*(0.5-sum2);//Initiation of 'e'
			
			}
      
				comp[array_index]=av_comp_alpha+_Complex_I*0.0;
				
		}
		fprintf(fp,"\n");		
	}
	for(INDEX=0; INDEX<max_particle; ++INDEX)
	{
	X = (double)(rand()/(RAND_MAX*1.0));
	Y = (double)(rand()/(RAND_MAX*1.0));
	
	X = nodes_x*X;
	Y = nodes_y*Y;
	R = r0;
	
	x[INDEX] = (int) X;
	y[INDEX] = (int) Y;
	rd[INDEX] = (int) R;
	iter=0;
		for (i=INDEX-1;i>=0;i--)
		{
		xxx=sqrt(pow((x[INDEX]-x[i]),2)+pow((y[INDEX]-y[i]),2));
		xxx1=sqrt(pow((x[INDEX]-x[i]+nodes_x),2)+pow((y[INDEX]-y[i]),2));
		xxx2=sqrt(pow((x[INDEX]-x[i]-nodes_x),2)+pow((y[INDEX]-y[i]),2));
		xxx3=sqrt(pow((x[INDEX]-x[i]),2)+pow((y[INDEX]-y[i]+nodes_y),2));
		xxx4=sqrt(pow((x[INDEX]-x[i]),2)+pow((y[INDEX]-y[i]-nodes_y),2));
		
			if(xxx<4*rd[INDEX] || xxx1<4*rd[INDEX] ||  xxx2<4*rd[INDEX]*ecc || xxx3<4*rd[INDEX]*ecc || xxx4<4*rd[INDEX]*ecc )
			//if(xxx<2.0*rd[INDEX])
			{
				
			X = (double)(rand()/(RAND_MAX*1.0));
			Y = (double)(rand()/(RAND_MAX*1.0));
	
			X = nodes_x*X;
			Y = nodes_x*Y;
			
	
			x[INDEX] = (int) X;
			y[INDEX] = (int) Y;
			
			i=INDEX;
			
			}
			
			continue;
		}
	for(i1=-nodes_x/2; i1< 3*nodes_x/2; ++i1)
		{
			for(i2=-nodes_y/2; i2< 3*nodes_y/2; ++i2)
			{
	 			if (i1<0)
				{
					i=i1+nodes_x;
				}
				else if (i1>=nodes_x)
				{
					i=i1-nodes_x;
				}
				else
				{
					i=i1;
				}
				if (i2<0)
				{
					j=i2+nodes_y;
				}
				else if (i2>=nodes_y)
				{
					j=i2-nodes_y;
				}
				else
				{
					j=i2;
				}
				tmp1=x[INDEX];
				tmp2=y[INDEX];
				if (pow((tmp1-i1)/(rd[INDEX]),2)+(pow((tmp2-i2)/(rd[INDEX]*ecc),2))<=1)
				{
					IJK=j+nodes_x*i;
					comp[IJK] =av_comp_beta+_Complex_I*0.0;
				
					for (num=0;num<=num_eta-1;num++)
					{
					eta1[num][array_index]=0.0+_Complex_I*0.0;
					}
				
				
				}
				
			}
		}
		//printf("%ld\n",iter);
	}
	for(j=0;j<nodes_y;j++)
	{
		for(i=0;i<nodes_x;i++)
		{
			array_index=i+nodes_x*j;
			
			//theta_c[array_index]=(comp[array_index]-av_comp_alpha)/(av_comp_beta-av_comp_alpha); //this is defined in a way such that it takes a value of 1 at the beta ppt. composition and a value of zero in tha matrix
			//theta_c[array_index]=(comp[array_index]-av_comp_alpha)/(av_comp_beta-av_comp_alpha); //this is defined in a way such that it takes a value of 1 at the beta ppt. composition and a value of zero in tha matrix
			theta_c[array_index]=comp[array_index];
			//theta_eta[array_index]=(eta[array_index])/av_eta; //this is defined in a way such that it takes a value of 1 at the beta ppt. eta and a value of zero in tha matrix
			theta_eta[array_index]=pow((eta1[0][array_index]),2.0);
			fi[array_index]=theta_eta[array_index]; //so fi is a function of eta

			fi_times_theta_eta[array_index]=fi[array_index]*theta_eta[array_index];

			fi_times_theta_c[array_index]=fi[array_index]*theta_c[array_index];
			comp_real=creal(comp[array_index]);
			eta_real=creal(eta[array_index]);
			
		}
		fprintf(fp,"\n");		
	}
	fclose(fp);
}
