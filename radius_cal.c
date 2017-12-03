//calculation of radius of the precipitate
double radius(fftw_complex *comp,fftw_complex *eta,int nodes_x,int nodes_y,double dx)

//#define dx 0.5
//#define Nx 512
//#define Ny 512
//int main ()
             {  
                int i,j,index,centre;
                FILE *fp,*fp1,*fp2;
                fftw_complex *comp_int,*eta_int,*comp_int;
                double x1,x2,radius; 
                  
                   // comp=(double*)malloc(Nx*Ny*sizeof(double));
	           // eta=(double*)malloc(Nx*Ny*sizeof(double));
                    comp_int=(fftw_complex*)malloc(Nx*sizeof(double));
	            eta_int=(fftw_complex*)malloc(Nx*sizeof(double));
                    

                //    fp = fopen("comp_and_eta_ep_c_0.0000ep_eta_0.0000_A_z_1.0_delta_1.0_supsat_0.00_timestep_20000.dat","r");

                  //  for(j=0;j<Ny;j++)
	          // {
		     // for (i=0;i<Nx;i++)
                                 //    {
                                   //     index=i*Nx+j;
                                     //   fscanf(fp,"%d%d%lf%lf",&i,&j,&comp[index],&eta[index]);
                                   //   }
                
	           //}
//fclose(fp);
//writing scanning data
               /*    fp1=fopen("scanned_data_300.dat","w");
       
        for (j=0;j<Ny;j++)
              {
                  for (i=0;i<Nx;i++)
                                  {
                                   index=i*Nx+j;
                                   fprintf(fp1,"%d\t%d\t%lf\t%lf\n",i,j,comp[index],eta[index]);
                                  }    
//fprintf(fp,"\n");          
              }
fclose(fp1);*/
//1 d profile generation
fp2=fopen("1d_20000_5.dat","w");
 //extracting the relevant profiles for the radius of ppt computation
	for(i=0;i<nodes_x;i++)
	{
		comp_int[i]=comp[i+nodes_x*(nodes_y/2)];
		eta_int[i]=eta[i+nodes_x*(nodes_y/2)];
	}

     for (i=0;i<Nx;i++)
           {
                  fprintf(fp2,"%d\t%lf\t%lf\n",i,comp_int[i],eta_int[i]);
              }    

fclose(fp2);
   // centre=Nx/2;
    //diff_comp=comp_int[centre]-0.5;
  for  (i=0;i<Nx/2;i++)
                 { 
              // diff_comp=comp_int[i]-0.5;
             //  diff_eta=eta_int[i]-0.5;
               if   (comp_int[i]<0.7&&comp_int[i+1]>0.7)
                     {
                    printf("%lf\t%lf\n%d\t%d\n",comp_int[i],comp_int[i+1],i,i+1);
		    //x1 = (i+1)- ((comp_int[i+1]-0.65)/(comp_int[i+1]-comp_int[i]));
                      x1=(2*i+1)*0.5;
		     printf("%lf\n",x1);
	
               	break;
                   }
}
/* for  (i=Nx/2;i<Nx;i++)
                 { 
              // diff_comp=comp_int[i]-0.5;
             //  diff_eta=eta_int[i]-0.5;
               if   (comp_int[i]>0.65&&comp_int[i+1]<0.65)
                     {
                     printf("%lf\t%lf\n",comp_int[i],comp_int[i+1]);
		     x2 = (i+1)- ((comp_int[i+1]-0.65)/(comp_int[i+1]-comp_int[i]));
		     printf("%lf\n",x2);
	
               	break;
                   }
}*/
                 radius=(Nx*0.5-x1)*0.5;
                
          printf("radius=%lf",radius);
}

