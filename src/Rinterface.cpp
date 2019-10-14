#include <RcppEigen.h>
#include <random>

#include "dist_functions.hpp"
#include "hmm_functions.hpp"

// [[Rcpp::depends(RcppEigen)]]

//' @param r vector of distances associatd with different BEFs
//' @param n_j matrix of integers denoting the start and length of each observations associated BEF distances
//' @param d a 1D grid of positive real values over which the differing intensities are evaluated
//' @param L component truncation number
//' @param K intensity cluster truncation number
//' @param J number of rows in r matrix; number of groups
// [[Rcpp::export]]
Rcpp::List fhmm_fit(
		const Eigen::ArrayXd &r,
        const Eigen::ArrayXXi &n_j,
        const Eigen::ArrayXd &d,
        const int &L,
        const int &K,
        const int &J,
		const int &iter_max,
		const int &warm_up,
		const int &thin,
		const int &seed,
		const int &chain,
		const int &num_posterior_samples
		){


    std::mt19937 rng;
    rng = std::mt19937(seed);

#include "during_sampling.hpp"
#include "sampling_containers.hpp"

	Eigen::ArrayXd alpha = Eigen::ArrayXd::Ones(K);
	Eigen::ArrayXXd theta = Eigen::ArrayXXd::Ones(L,K);
	double prior_mu = 0.0;
	double prior_sd = 1.0;
	std::normal_distribution<double> norm_draw(0,1);
	mus = rnorm(L,K,prior_mu,prior_sd,rng);
	pi = rdirichlet(alpha,rng);
	w = rdirichlet(theta,rng);
	double sample_ix = 0;
	Eigen::ArrayXXd mu_sum(L,K);
	mu_sum = Eigen::ArrayXd::Zero(L,K);



	/*for(int iter_ix = 1 ; iter_ix <= iter_max ; iter_ix ++){

		// update cluster assignment
		probs = dnorm_pi(J,r,n_j,pi,w,mus);

	
	for(int i = 0 ; i < J;  i++){
		prob = probs.row(i);
		std::discrete_distribution<int> d(prob.data(),prob.data() + prob.size());
		cluster_assignment(i) = d(rng);
	}

	component_probs = dnorm_w(r,n_j,w,mus,cluster_assignment);

		mu_sum = Eigen::ArrayXd::Zero(K);
		for(int k = 0 ; k < K; k ++){
			for(int i = 0; i < r.size(); i ++){
				mu_sum(k) += (cluster_assignment(i) == k)  ? r(i) : 0.0;
			}
		}

		// update mus
		for(int k = 0; k < K; k ++ ){
			cluster_count(k) = (cluster_assignment == k).count();
			mus(k) = norm_draw(rng)* sqrt( 1.0 / ( 1 + cluster_count(k) )  ) + mu_sum(k)  / ( 1 + cluster_count(k) ) ;
		}

		// update pis
		alpha_prime = alpha + cluster_count;
		pi = rdirichlet(alpha_prime,rng);

		if(iter_ix > warm_up && (iter_ix % thin == 0)){
			sampled_pis.row(sample_ix) = pi;
			sampled_mus.row(sample_ix) = mus;
			sample_ix += 1;
			for(int d_ = 0; d_ < d.size() ; d_ ++)
				density_samples(sample_ix,d_) = (pi * pow(2*PI,-.5) * exp(-.5 * pow(d(d_) - mus,2) )).sum() ;
		}
	}*/


    return(Rcpp::List::create(
				Rcpp::Named("density") = density_samples,
				Rcpp::Named("pi") = sampled_pis,
				Rcpp::Named("mu") = sampled_mus,
				Rcpp::Named("w") = w,
				Rcpp::Named("Cluster_Assignment") = cluster_assignment
				));
}
