
#ifndef  MIXTUREMODEL
#define  MIXTUREMODEL

class MixtureModel{

	public: 
		Eigen::ArrayXd r;

		MixtureModel(Eigen::ArrayXd &input_r){
		
			r = input_r;
		}

};

#endif 
