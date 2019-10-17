Eigen::ArrayXXd dnorm_pi(const int &J, const Eigen::ArrayXd &r, 
						 const Eigen::ArrayXXi &n_j, Eigen::ArrayXd &pi,
					     Eigen::ArrayXXd &w, Eigen::ArrayXXd &mu, 
               Eigen::ArrayXXd &taus){


	const int L = mu.rows();
	const int K = mu.cols();
	Eigen::ArrayXXd q(J,K);
	Eigen::ArrayXd tmp;

	for(int j = 0; j < J; j++){
		for(int k = 0; k < K; k ++){
			tmp = Eigen::ArrayXd::Zero(n_j(j,1));
			for(int l = 0; l < L; l ++){
			tmp = tmp + w(l,k) * pow(2 * M_PI * taus(l,k), -.5) * exp(- 1.0 / (2.0 * taus(l,k) ) * pow((r.segment(n_j(j,0),n_j(j,1)) - mu(l,k)),2 ));
			}
			q(j,k) = (pi(k)) * ((tmp)).prod();
		}
	}

  return(q);

}

Eigen::ArrayXXd dnorm_w(const Eigen::ArrayXd &r, const Eigen::ArrayXXi &n_j, 
                        Eigen::ArrayXXd &w,Eigen::ArrayXXd &mu,
                        Eigen::ArrayXXd &tau, Eigen::ArrayXi &zeta){

	const int L = mu.rows();
	const int J = zeta.size();
	Eigen::ArrayXXd b(r.size(),L);
  double log_factor = log(pow(10,-16)) + log(r.size());

	for(int l = 0; l < L; l ++){
		for(int j = 0; j< J; j++){
			for(int i = 0; i < n_j(j,1); i++){
				b(n_j(j,0)+i,l) = log(w(l,zeta(j)))  + R::dnorm(r(n_j(j,0) +i ), mu(l,zeta(j)), sqrt(tau(l,zeta(j))) , true);
        b(n_j(j,0)+i,l) = b(n_j(j,0)+i,l) > log_factor ? exp(b(n_j(j,0)+i,l)) : 0.0;
      }
		}
	}

	return(b);

}

Eigen::ArrayXXd dnorm(const Eigen::ArrayXd &r, Eigen::ArrayXd &pi, Eigen::ArrayXd &mu){

	Eigen::ArrayXXd probs(r.size(),pi.size());
	double log_factor = log(pow(10,-16)) + log(r.size());

	for(int k = 0; k < pi.size(); k ++){
		probs.col(k) = log(pi(k))   -.5 * pow((r-mu(k)),2);
		probs.col(k) = (probs.col(k) > log_factor).cast<double>() * exp(probs.col(k));
	}
	return(probs);
}
