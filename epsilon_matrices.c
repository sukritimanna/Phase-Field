//this is a routine for computing the epsilon matrices
void epsilon_mat(double epsilon_c[3][3],double epsilon_eta[3][3],double kron_delta[3][3],double epsilon_c_mag,double epsilon_eta_mag)
{
	int i,j;
	//creating the epsilon_c, epsilon_eta and kron_delta matrices
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			//assuming dilatational nature of the eigen strains
			if(i==j)
			{
				if(i!=2)
				{
					epsilon_c[i][j]=epsilon_c_mag;
					epsilon_eta[i][j]=epsilon_eta_mag;
				}
				else
				{
					epsilon_c[i][j]=0.0;
					epsilon_eta[i][j]=0.0;
				}
				kron_delta[i][j]=1.0;
			}
			else
			{
				epsilon_c[i][j]=0.0;
				epsilon_eta[i][j]=0.0;
				kron_delta[i][j]=0.0;
			}
		}		
	}	
	
}
