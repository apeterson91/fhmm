#ifndef DURING_SAMPLING
#define DURING SAMPLING


Eigen::ArrayXi cluster_assignment(J);
Eigen::ArrayXXd probs(J,K);
Eigen::ArrayXXd component_probs(r.size(),L);
Eigen::ArrayXd prob(K);
Eigen::ArrayXd pi(K);
Eigen::ArrayXXd w(L,K);
Eigen::ArrayXXd mus(L,K);
Eigen::ArrayXd cluster_count(K);
Eigen::ArrayXXd component_count(L,K);
Eigen::ArrayXd alpha_prime(K);
Eigen::ArrayXXd theta_prime(L,K);


#endif
