//this is a code which computes the homogeneous strain in the system 
void hom_str(double hom_strain[3][3],double average_periodic_stress[3][3],double average_eigen_stress[3][3],double kron_delta[3][3],double lambda_bar,double mu_bar,double mu_prime_bar,int iter)
{
	double trace_eigen_stress,trace_periodic_stress;	
	int i,j;
	//forming the traces of the tensors
	trace_eigen_stress=average_eigen_stress[0][0]+average_eigen_stress[1][1]+average_eigen_stress[2][2];
	trace_periodic_stress=average_periodic_stress[0][0]+average_periodic_stress[1][1]+average_periodic_stress[2][2];
	
	for(i=0;i<=2;i++)
	{
		for(j=0;j<=2;j++)
		{
			hom_strain[i][j]=lambda_bar*(trace_eigen_stress-trace_periodic_stress)*kron_delta[i][j]+2.0*mu_bar*(average_eigen_stress[i][j]-average_periodic_stress[i][j])+mu_prime_bar*((average_eigen_stress[0][0]-average_periodic_stress[0][0])*kron_delta[i][0]*kron_delta[j][0]+(average_eigen_stress[1][1]-average_periodic_stress[1][1])*kron_delta[i][1]*kron_delta[j][1]+(average_eigen_stress[2][2]-average_periodic_stress[2][2])*kron_delta[i][2]*kron_delta[j][2]);
			
		}			
	
	}
}
			
