#ifndef DURING_SAMPLING
#define DURING SAMPLING


Eigen::ArrayXi cluster_assignment(r.size());
Eigen::ArrayXXd probs(r.size(),K);
Eigen::ArrayXd prob(K);
Eigen::ArrayXd pi(K);
Eigen::ArrayXd mus(K);
Eigen::ArrayXd cluster_count(K);
Eigen::ArrayXd alpha_prime(K);


#endif
