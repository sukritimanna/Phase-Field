// this is a code for calculating the compliance tensor

void compl(double *S_11,double *S_12,double *S_44,double *lambda_bar,double *mu_bar, double *mu_prime_bar,double avg_C_11,double avg_C_12,double avg_C_44)

{
	*S_44=1.0/avg_C_44;
	*S_11=-(avg_C_11+avg_C_12)/((avg_C_11+2.0*avg_C_12)*(avg_C_12-avg_C_11));
	*S_12=avg_C_12/((avg_C_11+2.0*avg_C_12)*(avg_C_12-avg_C_11));
	//printf("The compliance components S_11=%lf\tS_12=%lf\tS_44=%lf\n",*S_11,*S_12,*S_44);
	*lambda_bar=*S_12;
	*mu_bar=*S_44;
	*mu_prime_bar=*S_11-*S_12-2.0*(*S_44);	
	//printf("The compliance tensor values we are going to use:lambda_bar=%lf\tmu_bar=%lf\tmu_prime_bar=%lf\n",*lambda_bar,*mu_bar,*mu_prime_bar);
}
