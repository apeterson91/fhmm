

#ifndef DIST_FUNCTIONS
#define DIST_FUNCTIONS

Eigen::ArrayXd rnorm( const int &n, double &mu, double &sd, std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);
	Eigen::ArrayXd out(n);

	for(int i = 0; i < n; i ++)
		out(i) = rnorm(rng)*sd + mu;

	return(out);
}


double rbeta(double &shape_1, double &shape_2, std::mt19937 &rng){

	std::chi_squared_distribution<double> rchisq_1(2 * shape_1);
	std::chi_squared_distribution<double> rchisq_2(2 * shape_2);

	double temp = rchisq_1(rng);
	double out = temp / (rchisq_2(rng) + temp);

	return(out);

}

Eigen::ArrayXd rbeta(int &n, double &shape_1, double &shape_2, std::mt19937 &rng){

	Eigen::ArrayXd out(n);
	for(int i = 0; i < n; i++)
		out(n) = rbeta(shape_1,shape_2,rng);

	return(out);
}



Eigen::ArrayXd rdirichlet(Eigen::ArrayXd &alpha, std::mt19937 &rng){


	Eigen::ArrayXd out(alpha.rows());

	for(int i = 0; i < alpha.rows() ; i++ ){
		std::gamma_distribution<double> rgamma(alpha(i),alpha(i));
		out(i) = rgamma(rng);
	}

	out = out / out.sum();
	return(out);
}


Eigen::ArrayXXd rdirichlet( int &n, Eigen::ArrayXd &alpha, std::mt19937  &rng){

	Eigen::ArrayXXd out(n,alpha.rows());

	for(int i = 0; i < n;  i ++)
		out.row(i) = rdirichlet(alpha,rng);

	return(out);
}

#endif
