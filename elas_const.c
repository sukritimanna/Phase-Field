//this routine is written for the purpose creating the elastic constants
void el_const(double C_44_mat,double delta,double nu,double A_z,double *C_12_mat, double *C_11_mat,double *C_12_ppt, double *C_11_ppt,double *C_44_ppt,double *lambda_mat,double *mu_mat,double *mu_prime_mat,double *lambda_ppt,double *mu_ppt,double *mu_prime_ppt,double *del_lambda,double *del_mu,double *del_mu_prime,double *del_C_12, double *del_C_11,double *del_C_44)
{

	//computing the elastic constants for the matrix
	(*C_12_mat)=(2.0*nu*C_44_mat)/(1-2.0*nu);
	(*C_11_mat)=(*C_12_mat)+((2.0*C_44_mat)/A_z);
	(*lambda_mat)=(*C_12_mat);
	(*mu_mat)=C_44_mat;
	(*mu_prime_mat)=(*C_11_mat)-(*C_12_mat)-2.0*C_44_mat;
	
	//computing the elastic constants for the precipitate
	(*C_44_ppt)=delta*C_44_mat;
	(*C_12_ppt)=(2.0*nu*(*C_44_ppt))/(1-2.0*nu);
	(*C_11_ppt)=(*C_12_ppt)+((2.0*(*C_44_ppt))/A_z);
	(*lambda_ppt)=(*C_12_ppt);
	(*mu_ppt)=(*C_44_ppt);
	(*mu_prime_ppt)=(*C_11_ppt)-(*C_12_ppt)-2.0*(*C_44_ppt);
	

	//computing the del_lambda,del_mu and del_mu_prime
	(*del_lambda)=(*lambda_ppt)-(*lambda_mat);
	(*del_mu)=(*mu_ppt)-(*mu_mat);
	(*del_mu_prime)=(*mu_prime_ppt)-(*mu_prime_mat);
	
	
	//computing del_C_11, del_C_12 and del_C_44
	(*del_C_11)=(*C_11_ppt)-(*C_11_mat);
	(*del_C_12)=(*C_12_ppt)-(*C_12_mat);
	(*del_C_44)=(*C_44_ppt)-(C_44_mat); 
}

