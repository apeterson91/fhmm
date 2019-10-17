#include <RcppEigen.h>
#include <random>

#include "dist_functions.hpp"
#include "hmm_functions.hpp"
#include "auxiliary_functions.hpp"

// [[Rcpp::depends(RcppEigen)]]

//' @param r vector of distances associatd with different BEFs
//' @param n_j matrix of integers denoting the start and length of each observations associated BEF distances
//' @param d a 1D grid of positive real values over which the differing intensities are evaluated
//' @param L number of cluster mixture components
//' @param K number of function components
//' @param J number of rows in r matrix; number of groups
//' @param mu_0 prior mean for mean parameter
//' @param kappa_0 prior variance parameter for mean normal prior
//' @param nu_0 prior degrees of freedom for variance Inv Chisq prior
//' @param sigma_0 prior scale for variance Inv-Chisq distribtuion
//' @param iter_max number of total iterations to run the sampler for
//' @param warm_up number of iterations to discard as burn-in or warm_up
//' @param thin number of iterations to thin posterior sample draws by
//' @param seed seed with which to initialize random number generator
//' @param chain used for labeling
//' @param num_posterior_samples number of posterior samples kept after thinning
// [[Rcpp::export]]
Rcpp::List fhmm_fit(
        const Eigen::ArrayXd &r,
        const Eigen::ArrayXXi &n_j,
        const Eigen::ArrayXd &d,
        const int &L,
        const int &K,
        const int &J,
        const double &mu_0,
        const double &kappa_0,
        const int &nu_0,
        const double &sigma_0,
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
	int d_length = d.size();
	std::normal_distribution<double> norm_draw(0,1);
  taus = initialize_tau(L,K,nu_0,sigma_0);
	mus = rnorm(L,K,prior_mu,prior_sd,rng);
  mus = mus * sqrt(taus); 
	pi = rdirichlet(alpha,rng);
	w = rdirichlet(theta,rng);
	double sample_ix = 0;
  double s_n = 0;
  double mu_n = 0;
	Eigen::ArrayXXd mu_sum(L,K);
  Eigen::ArrayXXd tau_sum(L,K);



	for(int iter_ix = 1 ; iter_ix <= iter_max ; iter_ix ++){
    print_progress(iter_ix, warm_up, iter_max, chain);

		// update cluster assignment
		probs = dnorm_pi(J,r,n_j,pi,w,mus,taus);

		for(int j = 0 ; j < J; j ++){
			prob = probs.row(j);
			std::discrete_distribution<int> d(prob.data(),prob.data() + prob.size());
			cluster_assignment(j) = d(rng);
		}

		component_probs = dnorm_w(r,n_j,w,mus,taus,cluster_assignment);

		component_count = Eigen::ArrayXXd::Zero(L,K);
		mu_sum = Eigen::ArrayXXd::Zero(L,K);
    tau_sum = Eigen::ArrayXXd::Zero(L,K);

		//assign distances to components within clusters
		for(int j = 0; j < J; j++){
			for(int i = 0; i < n_j(j,1) ; i ++){
				for(int l = 0; l < L; l++){
					c_prob(l) = component_probs(n_j(j,0)+i,l);
				}
				std::discrete_distribution<int> d(c_prob.data(),c_prob.data() + c_prob.size());
				component_assignment(n_j(j,0)+i) = d(rng);
				component_count(component_assignment(n_j(j,0)+i),cluster_assignment(j)) += 1;
				mu_sum(component_assignment(n_j(j,0)+i),cluster_assignment(j)) += r(n_j(j,0)+i);
        tau_sum(component_assignment(n_j(j,0)+i),cluster_assignment(j)) += pow(r(n_j(j,0)+i),2);
			}
		}

		for(int k = 0; k < K; k ++ )
				cluster_count(k) = (cluster_assignment == k).count();

		// update mus
		for(int k = 0; k < K; k ++ ){
				cluster_count(k) = (cluster_assignment == k).count();
			for(int l = 0; l < L; l ++){
        if(component_count(l,k)==0){
          mus(l,k) = norm_draw(rng)*prior_sd + prior_mu;
          taus(l,k) = nu_0 * sigma_0 / R::rchisq(1);
        }
        else{
          s_n  = nu_0 * sigma_0 + (tau_sum(l,k) - (pow(mu_sum(l,k),2) / component_count(l,k) ) ) + (kappa_0 * component_count(l,k) / (kappa_0 + component_count(l,k))) * pow( (mu_sum(l,k) / component_count(l,k) - mu_0 ),2);
          taus(l,k) = s_n / R::rchisq(nu_0 + component_count(l,k));
          mu_n = ( (kappa_0  / taus(l,k)) * mu_0   + mu_sum(l,k) / taus(l,k)  ) / ( kappa_0 / taus(l,k) + component_count(l,k) / taus(l,k));
          s_n = 1.0 /( (kappa_0 /taus(l,k)) + (component_count(l,k) / taus(l,k)) ) ;
          mus(l,k) =  norm_draw(rng) * sqrt(s_n) + mu_n;
        }
      }
		}

		// update pis, ws
		alpha_prime = alpha + cluster_count;
		theta_prime = theta + component_count;
		pi = rdirichlet(alpha_prime,rng);
		w = rdirichlet(theta_prime,rng);

    // store samples
		if(iter_ix > warm_up && (iter_ix % thin == 0)){
			cluster_assignments.row(sample_ix) = cluster_assignment;
			Eigen::Map<Eigen::RowVectorXd> mu_samp(mus.data(),mus.size());
			Eigen::Map<Eigen::RowVectorXd> tau_samp(taus.data(),taus.size());
			Eigen::Map<Eigen::RowVectorXd> w_samp(w.data(),w.size());
			sampled_ws.row(sample_ix) = w_samp;
			sampled_mus.row(sample_ix) = mu_samp;
      sampled_taus.row(sample_ix) = tau_samp;
			sampled_pis.row(sample_ix) = pi;
			for(int k = 0; k < K; k++){
				for(int d_ix = 0; d_ix < d_length; d_ix ++){
					for(int l = 0; l < L; l++){
						intensities(sample_ix, k*d_length + d_ix) += w(l,k) * R::dnorm(d(d_ix),mus(l,k),1,false);
           }
				 global_intensity(sample_ix,d_ix) += pi(k) *  intensities(sample_ix,k*d_length+d_ix);
       }
      }
			sample_ix += 1;
		}

	}

    return(Rcpp::List::create(
				Rcpp::Named("pi") = sampled_pis,
				Rcpp::Named("mu") = sampled_mus,
        Rcpp::Named("tau") = sampled_taus,
        Rcpp::Named("w") = sampled_ws,
				Rcpp::Named("Cluster_Assignment") = cluster_assignments,
				Rcpp::Named("global_intensity") = global_intensity,
				Rcpp::Named("cluster_intensity") = intensities
				));
}
