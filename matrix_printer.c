//this is a code which will be used for printing the stress, strain or any such 3X3 matrices
void prt(double mat[3][3])
{
	int i,j;
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			printf("%lf\t",mat[i][j]);
		}
		printf("\n");
	}
}
