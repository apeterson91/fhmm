
#ifndef SAMPLING_CONTAINERS
#define SAMPLING_CONTAINERS

Eigen::ArrayXXi cluster_assignments(num_posterior_samples,J);
Eigen::MatrixXd sampled_mus(num_posterior_samples,L * K);
Eigen::MatrixXd sampled_taus(num_posterior_samples,L * K);
Eigen::MatrixXd sampled_ws(num_posterior_samples,L*K);
Eigen::ArrayXXd sampled_pis(num_posterior_samples,K);
Eigen::ArrayXXd intensities(num_posterior_samples,K * d.size());
Eigen::ArrayXXd global_intensity(num_posterior_samples,d.size());

cluster_assignments = Eigen::ArrayXXi::Zero(num_posterior_samples,J);
sampled_mus = Eigen::ArrayXXd::Zero(num_posterior_samples, L * K);
sampled_taus = Eigen::ArrayXXd::Zero(num_posterior_samples, L * K);
sampled_ws = Eigen::ArrayXXd::Zero(num_posterior_samples,L * K);
sampled_pis = Eigen::ArrayXXd::Zero(num_posterior_samples,K);
intensities = Eigen::ArrayXXd::Zero(num_posterior_samples,K * d.size());
global_intensity = Eigen::ArrayXXd::Zero(num_posterior_samples, d.size());

#endif
