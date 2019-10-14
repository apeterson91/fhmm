
#ifndef SAMPLING_CONTAINERS
#define SAMPLING_CONTAINERS

Eigen::ArrayXXd cluster_assignments(num_posterior_samples,r.size());
Eigen::ArrayXXd sampled_mus(num_posterior_samples,L*K);
Eigen::ArrayXXd sampled_pis(num_posterior_samples,K);
Eigen::ArrayXXd density_samples(num_posterior_samples,d.size());


#endif
