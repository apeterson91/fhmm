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
        const Eigen::MatrixXi &n_j,
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
	double prior_mu = 0.0;
	double prior_sd = 1.0;
	std::normal_distribution<double> norm_draw(0,1);
	mus = rnorm(K,prior_mu,prior_sd,rng);
	pi = rdirichlet(alpha,rng);
	

	for(int iter_ix = 1 ; iter_ix <= iter_max ; iter_ix ++){

		// update cluster assignment
		probs = dnorm(r,pi,mus);
		for(int i = 0 ; i < r.size(); i++){
			prob = probs.row(i);
			std::discrete_distribution<int> d(prob.data(),prob.data() + prob.size());
			cluster_assignment(i) = d(rng);
		}

		for(int k = 0; k < K ; k ++){
			cluster_count(k) = (cluster_assignment == k).count();
		}

		// update mus
		for(int k = 0; k < K; k ++ ){
			mus(k) = norm_draw(rng) + ((cluster_assignment == k).cast<double>() * r).sum()  / ( 1 + cluster_count(k) ) ;
		}
		// update pis
		alpha_prime = alpha + cluster_count;
		pi = rdirichlet(alpha_prime,rng);

	
	}


    return(Rcpp::List::create(Rcpp::Named("cluster_assignment") = pi));
}
