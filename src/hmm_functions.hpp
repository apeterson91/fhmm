
Eigen::ArrayXXd dnorm(const Eigen::ArrayXd &r, Eigen::ArrayXd &pi, Eigen::ArrayXd &mu){

	Eigen::ArrayXd probs(r.size(),pi.size());

	for(int k = 0; k < pi.size(); k ++){
		probs.col(k) = pi(k) * exp( -.5 * (r - mu(k)).matrix().dot((r-mu(k)).matrix()) );
	}

	return(probs);
}
