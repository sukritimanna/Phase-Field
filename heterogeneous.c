//this is a routine to calculate the heterogeneous strains, again we have used general names for the displacement arrays passed  

void het_strains(fftw_complex *ux,fftw_complex *uy,fftw_complex *het11,fftw_complex *het22,fftw_complex *het12,fftw_plan plb,double *kx,double *ky,int nodes_x,int nodes_y,int iter)
{
	//FILE *fp;
	//char fname[300];
	//sprintf(fname,"../DATA/het_first_inhomogen");
	//if((fp=fopen(fname,"w"))==NULL)
	//{
	//	printf("het homogeneous can't be written\n");
	//	exit(1);
	//}
	int i,j;
	double sum11=0.0,sum12=0.0,sum22=0.0;
	
	//populating the fourier transformed values of the periodic strain
	//here we have enforced our understanding of the dimensions as indices in tensors
	for(j=0;j<nodes_y;j++)
	{
		for(i=0;i<nodes_x;i++)
		{
			het11[i+nodes_x*j]=kx[i]*_Complex_I*ux[i+nodes_x*j];
			het22[i+nodes_x*j]=ky[j]*_Complex_I*uy[i+nodes_x*j];
			//this definition is in accordance with the one used in the formulation
			het12[i+nodes_x*j]=0.5*(kx[i]*_Complex_I*uy[i+nodes_x*j]+ky[j]*_Complex_I*ux[i+nodes_x*j]); 
		}
	}
	//taking the inverse fourier transform of the arrays
	fftw_execute_dft(plb,het11,het11);
	fftw_execute_dft(plb,het22,het22);
	fftw_execute_dft(plb,het12,het12);
	//normalising the arrays
	for(j=0;j<nodes_y;j++)
	{
		for(i=0;i<nodes_x;i++)
		{
			het11[i+nodes_x*j]=het11[i+nodes_x*j]/(double)(nodes_x*nodes_y);
			het22[i+nodes_x*j]=het22[i+nodes_x*j]/(double)(nodes_x*nodes_y);
			het12[i+nodes_x*j]=het12[i+nodes_x*j]/(double)(nodes_x*nodes_y);
			//fprintf(fp,"%d\t%d\t%lf\t%lf\t%lf\n",j,i,creal(het11[i+nodes_x*j]),creal(het22[i+nodes_x*j]),creal(het12[i+nodes_x*j]));
			sum11=sum11+creal(het11[i+nodes_x*j]);
			sum22=sum22+creal(het22[i+nodes_x*j]);
			sum12=sum12+creal(het12[i+nodes_x*j]);
		}
	}
	//fclose(fp);
	if(fabs(sum11)>pow(10,-8) || fabs(sum22)>pow(10,-8) || fabs(sum12)>pow(10,-8))
	{
		printf("The summations of heterogeneous strains over the entire domain have a problem\n");
		printf("the three sums are sum11=%lf,sum22=%lf,sum12=%lf\n",sum11,sum22,sum12);
		printf("this happens at iteration no. %d\n",iter);
		exit(1);
	}
}
