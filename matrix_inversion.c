//doing a 3x3 matrix inversion by the gauss jordan elimination technique
//though the variables in this routine reflect its application in inverting omega matrices it is also used to invert the contracted modulus matrices 
//the expressions used in this routine are ok 
//we are going to tailor this routine for inverting a 3X3 matrix
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void matrix_invert(double omg[3][3],double omg_in[3][3])
{
	double eliminate[3][3],sum1,sum2,pro1[3][3],pro2[3][3],replace[3];
	 
	int i,j,k,t,s;

	//initialising the inverted matrix as an identity matrix
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			if(i==j)
				omg_in[i][j]=1.0;
			else
				omg_in[i][j]=0.0;
				
		}
	}
	//as the matrix is NxN we need N+1 iterations to arrive at the inverse
	//the first N iterations are for removing the non-diagonal elements and the last one is for scaling the diagonal elements
	for(t=0;t<4;t++)
	{
		if(t==3)
		{
			//the last iteration
			//creating the elimination matrix
			for(i=0;i<3;i++)
			{
				for(j=0;j<3;j++)
				{
					if(i==j)
						eliminate[i][j]=1.0/omg[i][j];
					else
						eliminate[i][j]=0.0;
				}
			}
		}
		else
		{
			s=t;
			//partial pivoting done if necessary
			while(s<=1)
			{
				s++;
				if(fabs(omg[t][t])<fabs(omg[s][t]))
				{
					for(i=0;i<3;i++)
					{	
						replace[i]=omg[t][i];
						omg[t][i]=omg[s][i];
						omg[s][i]=replace[i];
						replace[i]=omg_in[t][i];
						omg_in[t][i]=omg_in[s][i];
						omg_in[s][i]=replace[i];
					}
				}
			}
			//populating the elimination matrix
			for(i=0;i<3;i++)
			{
				for(j=0;j<3;j++)
				{	
					if(i==j)
						eliminate[i][j]=1.0;
					else if((j==t) && (i!=j))
						eliminate[i][j]=-(omg[i][j]/omg[j][j]);
					else
						eliminate[i][j]=0.0;
				}
			}
		}
		//multiplying the elimination matrix with the original matrix and with the identity matrix
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				sum1=0.0;
				sum2=0.0;
				for(k=0;k<3;k++)
				{
					sum1=sum1+(eliminate[i][k]*omg[k][j]);
					sum2=sum2+(eliminate[i][k]*omg_in[k][j]);
				}
				pro1[i][j]=sum1;
				pro2[i][j]=sum2;
			}
		}
		//creating the new matrices
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				omg[i][j]=pro1[i][j];
				omg_in[i][j]=pro2[i][j];
			}
		}
	}
	
}
