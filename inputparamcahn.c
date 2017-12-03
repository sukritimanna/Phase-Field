//this is a routine to read the system parameters

void inputparam(double *sys_size_x,double *sys_size_y,double *dx,double *dy,double *M,double *L,double *K_comp,double *K_eta,double *epsilon_c_mag,double *epsilon_eta_mag,double *C_44_mat,double *nu,double *A_z,double *delta,double *radius,double *av_comp_alpha,double *av_comp_beta,double *av_eta,double *toler,double *dt,int *total_time,int *output_time_interval,int *num_eta,double *A, double *B, double *Z)
{
	char str[100],fname[1000];
	FILE *fpread,*fpwrite;
	//opening the files for the reading and writing of system parameters
	if((fpread=fopen("../INPAR/read.dat","r"))==NULL)
	{
		printf("System parameters can't be read\n");
		exit(1);
	}
	//reading from one and writing in the other
	fscanf(fpread,"%s%lf",str,sys_size_x);
	fscanf(fpread,"%s%lf",str,sys_size_y);
	fscanf(fpread,"%s%lf",str,dx);
	fscanf(fpread,"%s%lf",str,dy);
	fscanf(fpread,"%s%lf",str,M);
	fscanf(fpread,"%s%lf",str,L);
	fscanf(fpread,"%s%lf",str,K_comp);
	fscanf(fpread,"%s%lf",str,K_eta);
	fscanf(fpread,"%s%lf",str,epsilon_c_mag);
	fscanf(fpread,"%s%lf",str,epsilon_eta_mag);
	fscanf(fpread,"%s%lf",str,C_44_mat);
	fscanf(fpread,"%s%lf",str,nu);
	fscanf(fpread,"%s%lf",str,A_z);
	fscanf(fpread,"%s%lf",str,delta);
	fscanf(fpread,"%s%lf",str,radius);
	fscanf(fpread,"%s%lf",str,av_comp_alpha);
	fscanf(fpread,"%s%lf",str,av_comp_beta);
	fscanf(fpread,"%s%lf",str,av_eta);
	fscanf(fpread,"%s%lf",str,toler);
	fscanf(fpread,"%s%lf",str,dt);
	fscanf(fpread,"%s%d",str,total_time);
	fscanf(fpread,"%s%d",str,output_time_interval);
	fscanf(fpread,"%s%d",str,num_eta);
	fscanf(fpread,"%s%lf",str,A);
	fscanf(fpread,"%s%lf",str,B);
	fscanf(fpread,"%s%lf",str,Z);
	//fscanf(fpread,"%s%lf",str,sup_sat);
	//fscanf(fpread,"%s%lf",str,temp);
	//fscanf(fpread,"%s%lf",str,delta_f);
	sprintf(fname,"../DATA/param_ep_c_%.4lfep_eta_%.4lf_A_z_%.1lf_delta_%.1lf.dat",*epsilon_c_mag,*epsilon_eta_mag,*A_z,*delta);
	if((fpwrite=fopen(fname,"w"))==NULL)
	{
		printf("System parameters can't be written\n");
		exit(1);
	}
	fprintf(fpwrite,"sys_size_x\t%lf\n",*sys_size_x);
	fprintf(fpwrite,"sys_sixe_y\t%lf\n",*sys_size_y);
	fprintf(fpwrite,"dx\t%lf\n",*dx);
	fprintf(fpwrite,"dy\t%lf\n",*dy);
	fprintf(fpwrite,"M\t%lf\n",*M);
	fprintf(fpwrite,"L\t%lf\n",*L);
	fprintf(fpwrite,"K_comp\t%lf\n",*K_comp);
	fprintf(fpwrite,"K_eta\t%lf\n",*K_eta);
	fprintf(fpwrite,"epsilon_c_mag\t%lf\n",*epsilon_c_mag);
	fprintf(fpwrite,"epsilon_eta_mag\t%lf\n",*epsilon_eta_mag);
	fprintf(fpwrite,"C_44_mat\t%lf\n",*C_44_mat);
	fprintf(fpwrite,"nu\t%lf\n",*nu);
	fprintf(fpwrite,"A_z\t%lf\n",*A_z);
	fprintf(fpwrite,"delta\t%lf\n",*delta);
	fprintf(fpwrite,"radius\t%lf\n",*radius);
	fprintf(fpwrite,"av_comp_alpha\t%lf\n",*av_comp_alpha);
	fprintf(fpwrite,"av_comp_beta\t%lf\n",*av_comp_beta);
	fprintf(fpwrite,"av_eta\t%lf\n",*av_eta);
	fprintf(fpwrite,"toler\t%.12lf\n",*toler);
	fprintf(fpwrite,"dt\t%.12lf\n",*dt);
	fprintf(fpwrite,"total_time\t%d\n",*total_time);
	fprintf(fpwrite,"output_time_interval\t%d\n",*output_time_interval);
	fprintf(fpwrite,"num_eta\t%d\n",*num_eta);
	fprintf(fpwrite,"A\t%lf\n",*A);
	fprintf(fpwrite,"B\t%lf\n",*B);
	fprintf(fpwrite,"Z\t%lf\n",*Z);
	//fprintf(fpwrite,"sup_sat\t%lf\n",*sup_sat);
	//fprintf(fpwrite,"temp\t%lf\n",*temp);
	//fprintf(fpwrite,"delta_f\t%lf\n",*delta_f);
	fclose(fpread);
	fclose(fpwrite);
}				
