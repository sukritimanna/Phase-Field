//this function computes the stress considering cubic anisotropy using the shortened formulation

void str(double sigma[3][3],double epsilon[3][3],double kron_delta[3][3],double lambda,double mu,double mu_prime)
{
	int i,j,k,l;
	double trace_epsilon;
	trace_epsilon=epsilon[0][0]+epsilon[1][1]+epsilon[2][2];
	for(i=0;i<=2;i++)
	{
		for(j=0;j<=2;j++)
		{
			sigma[i][j]=lambda*trace_epsilon*kron_delta[i][j]+2.0*mu*epsilon[i][j]+mu_prime*(epsilon[0][0]*kron_delta[i][0]*kron_delta[j][0]+epsilon[1][1]*kron_delta[i][1]*kron_delta[j][1]+epsilon[2][2]*kron_delta[i][2]*kron_delta[j][2]);
			
		}			
		
	}
}
