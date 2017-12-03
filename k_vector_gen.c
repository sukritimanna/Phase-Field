//this routine is written for the generation of k vectors

void k_vec(double *kx, double *ky, double *sq_mag_k_vectors, int nodes_x,int nodes_y,double sys_size_x,double sys_size_y)
{
	int half_mark_x,half_mark_y,i,j;
	long array_index;
//populating the array for the k-vectors
	half_mark_x=nodes_x/2;
	half_mark_y=nodes_y/2;
	for(j=0;j<nodes_y;j++)
	{
		//in the following calculations kx and ky have been expressed in the most appropriate way to negate cancellation
		if(j>half_mark_y)
			ky[j]=(2.0*M_PI*(double)(j-nodes_y))/sys_size_y;
		else
			ky[j]=(2.0*M_PI*(double)(j))/sys_size_y;
		for(i=0;i<nodes_x;i++)
		{
			array_index=i+nodes_x*j;
			if(i>half_mark_x)
			 	kx[i]=(2.0*M_PI*(double)(i-nodes_x))/sys_size_x;
			else
				kx[i]=(2.0*M_PI*(double)(i))/sys_size_x;
			sq_mag_k_vectors[array_index]=(kx[i]*kx[i]+ky[j]*ky[j]);		
		}
	} 
}
