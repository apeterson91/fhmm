
Eigen::ArrayXXd dnorm(const Eigen::ArrayXd &r, Eigen::ArrayXd &pi, Eigen::ArrayXd &mu){

	Eigen::ArrayXXd probs(r.size(),pi.size());
	double log_factor = log(pow(10,-16)) + log(r.size());

	for(int k = 0; k < pi.size(); k ++){
		probs.col(k) = log(pi(k))   -.5 * pow((r-mu(k)),2); 
		probs.col(k) = (probs.col(k) > log_factor).cast<double>() * exp(probs.col(k));
	}
	return(probs);
}
