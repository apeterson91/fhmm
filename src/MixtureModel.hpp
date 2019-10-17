
#ifndef  MIXTUREMODEL
#define  MIXTUREMODEL

class MixtureModel{

	public: 
		Eigen::ArrayXd r;
    Eigen::ArrayXXi n_j;
    Eigen::ArrayXXi prior_alpha;
    Eigen::ArrayXXi prior_theta;
    Eigen::ArrayXXi alpha;
    Eigen::ArrayXXi theta;
    Eigen::ArrayXXd mu
    int J;
    int L;
    int K;


		MixtureModel(const Eigen::ArrayXd &input_r,
                 const Eigen::ArrayXXi &input_n,
                 const int &input_J,
                 const int &input_L
                 const int &input_K,
                 const double &prior_mu,
                 const double &prior_sd,
                 const Eigen::ArrayXXd &prior_alpha,
                 const Eigen::ArrayXd &prior_theta
                 ){
		
			r = input_r;
      n_j = input_n; 
      J = input_J;
      K = input_K;
      L = input_L;
		}

};

#endif 
