//this is a code which computes the average of the elastic constants over the entire domain
void avg(double *avg_C_11,double *avg_C_12,double *avg_C_44,int nodes_x,int nodes_y,double dx,double dy,double C_44_mat,double del_C_44,double nu,double A_z,fftw_complex *fi)
{
	int i,j,count=0;
	double block_val;
	long array_index,neighbour_right,neighbour_up,neighbour_diagonal;
	*avg_C_44=0.0;
	//summming the C_44 from all the nodes (we are using the trapezoidal rule of numerical integration which means we go to a particular block average out the values from all the surrounding nodes and then multipliy that with the area of the block)
	for(j=0;j<nodes_y;j++)
	{
		for(i=0;i<nodes_x;i++)
		{
			//block_val=0.0;
			array_index=i+nodes_x*j;
			//neighbour_right=(i+1)+nodes_x*j;
			//neighbour_up=i+nodes_x*(j+1);
			//neighbour_diagonal=(i+1)+nodes_x*(j+1);
			(*avg_C_44)+=C_44_mat+del_C_44*creal(fi[array_index]);
			
			//*avg_C_44=(*avg_C_44)+block_val*dx*dy;
			//count++;
		}
	}
	//creating the average (by dividing the sum obtained in the above step by the total no. of blocks)
	*avg_C_44=(*avg_C_44)/(nodes_x*nodes_y);
	//if(count!=((nodes_x-1)*(nodes_y-1)))
	//{
	//	printf("Looping incorrect in average.c routine\n");
	//	exit(1);
	//}
	//printf("The average C_44 of the system is=%lf\n",*avg_C_44);
	//computing the avg_C_12 and avg_C_11 of the system by using the avg_C_44 of the system
	*avg_C_12=(2.0*nu*(*avg_C_44))/(1-2.0*nu);
	*avg_C_11=(*avg_C_12)+((2.0*(*avg_C_44))/A_z);
//	printf("The average C_11 of the system=%lf\n",*avg_C_11);	
//	printf("The average C_12 of the system=%lf\n",*avg_C_12);
}



