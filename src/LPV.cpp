/*
 * Copyright (C) 2016 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * Author: Sina Mirrazavi
 * email:   sina.mirrazavi@epfl.ch
 * website: lasa.epfl.ch
 *
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 */


#include "LPV.h"


bool is_file_exist(const char *fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}


const double PI = 3.14159265358979323846264338327950288419716939937510;

void LPV::initialize(int Num_Com,int Num_state)
{

	/* Declare the number of the components
	 * and the dimension  of the state  */

	Num_Com_=Num_Com;
	Num_State_=Num_state;

	A_Matrix_ =new MatrixXd[Num_Com_];	for(int s=0; s<Num_Com_; s++ ){A_Matrix_[s].resize(Num_State_,Num_State_);}
	Prior_ =new double[Num_Com_];
	Mu_ =new VectorXd[Num_Com_];		for(int s=0; s<Num_Com_; s++ ){	Mu_[s].resize(Num_State_);	}
	Sigma_ =new MatrixXd[Num_Com_];		for(int s=0; s<Num_Com_; s++ ){	Sigma_[s].resize(Num_State_,Num_State_);	}

	Theta_.resize(Num_Com_);			Theta_.setZero();

}
void LPV::initialize_A(const char  *path_){

	/* Initialize A
	 * path_ is the path of A matrix*/


	MatrixXd fMatrix(1,1);fMatrix.setZero();


	if (is_file_exist(path_))
	{
		fMatrix=readMatrix(path_);
	}
	else
	{
		cout<<"The provided path does not exist."<<endl;
		cout<<"path_A "<<endl;
		cout<<path_<<endl;
		ERROR();
	}


	if ((fMatrix.rows()!=Num_Com_*Num_State_)||(fMatrix.cols()!=Num_State_))
	{
		cout<<"Initialization of the A matrices is wrong."<<endl;
		cout<<"A_V: "<<endl;cout<<fMatrix<<endl;
		cout<<"Number of Components is: "<<Num_Com_<<endl;
		cout<<"Dimension of states is: "<<Num_State_<<endl;
		cout<<"Dimension of states of A_V is: "<<fMatrix.cols()<<endl;
		cout<<"Path is "<<path_<<endl;
		ERROR();
	}
	int j = 0;



	for(int s=0; s<Num_Com_; s++ ){
		for(int i=0; i<Num_State_; i++ ){
			A_Matrix_[s].row(i)=fMatrix.row(j);
			j++;
		}
	}
	for(int s=0; s<Num_Com_; s++ ){
		cout<<"A_Matrix["<<s<<"]"<<endl;
		cout<<A_Matrix_[s]<<endl;
	}

}
void LPV::initialize_theta(const char  *path_prior_,const char  *path_mu_,const char  *path_sigma_){

	/* Initialize scheduling parameters
	 *  path_prior_ is the path of the prior matrix
	 *	path_mu_ is the path of the mean matrix
	 *	path_sigma_ is the path of the covariance matrix */
	cout<<"Initializing the model of LPV"<<endl;
	MatrixXd fMatrix;
	if (is_file_exist(path_prior_))
	{
		fMatrix=readMatrix(path_prior_);
	}
	else
	{
		cout<<"The provided path does not exist."<<endl;
		cout<<"path_prior_LPV "<<endl;
		cout<<path_prior_<<endl;
		ERROR();
	}
	if ((fMatrix.cols()!=Num_Com_)||(fMatrix.rows()!=1))
	{
		cout<<"Initialization of Prior is wrong."<<endl;
		cout<<"Number of components is: "<<Num_Com_<<endl;
		cout<<"Dimension of states of Prior is: "<<fMatrix.cols()<<endl;
		ERROR();
	}

	cout<<"Prior"<<endl;
	for (int i=0; i<Num_Com_; i++)
	{
		Prior_[i]=fMatrix(0,i);;
	}

	for (int i=0; i<Num_Com_; i++)
	{
		cout<<Prior_[i]<<endl;
	}


	fMatrix.setZero();
	if (is_file_exist(path_mu_))
	{
		fMatrix=readMatrix(path_mu_);
	}
	else
	{
		cout<<"The provided path does not exist."<<endl;
		cout<<"path_mu_LPV "<<endl;
		cout<<path_mu_<<endl;
		ERROR();
	}


	if ((fMatrix.cols()!=Num_Com_)||(fMatrix.rows()!=Num_State_))
	{
		cout<<"Initialization of Mean is wrong."<<endl;
		cout<<"Number of components is: "<<Num_Com_<<endl;
		cout<<"Dimension of states is: "<<Num_State_<<endl;
		cout<<"Dimension of states of Mean is: "<<fMatrix.rows()<<"*"<<fMatrix.cols()<<endl;
		ERROR();
	}

	for(int s=0; s<Num_Com_; s++ )
	{
		Mu_[s]=fMatrix.col(s);

	}

	for(int s=0; s<Num_Com_; s++ )
	{
		cout<<"Mu["<<s<<"]."<<endl;
		cout<<Mu_[s]<<endl;

	}


	fMatrix.resize(1,1);fMatrix.setZero();

	if (is_file_exist(path_sigma_))
	{
		fMatrix=readMatrix(path_sigma_);
	}
	else
	{
		cout<<"The provided path does not exist."<<endl;
		cout<<"path_sigma_LPV "<<endl;
		cout<<path_sigma_<<endl;
		ERROR();
	}

	if ((fMatrix.rows()!=Num_Com_*Num_State_)||(fMatrix.cols()!=Num_State_))
	{
		cout<<"Initialization of the covariance matrix is wrong."<<endl;
		cout<<"the covariance matrix : "<<endl;cout<<fMatrix<<endl;
		cout<<"Number of components is: "<<Num_Com_<<endl;
		cout<<"Dimension of states is: "<<Num_State_<<endl;
		cout<<"Dimension of states of the covariance matrix is: "<<fMatrix.rows()<<"*"<<fMatrix.cols()<<endl;
		ERROR();
	}

	int s=0;

	for(int i=0; i<Num_Com_; i++ ){
		for(int j=0; j<Num_State_; j++ ){
			Sigma_[i].row(j)=fMatrix.row(s);
			s++;
		}
	}

	for(int i=0; i<Num_Com_; i++ ){
		cout<<"Sigma["<<i<<"]."<<endl;
		cout<<Sigma_[i]<<endl;
	}

}
MatrixXd LPV::Calculate_A(VectorXd X){

	/* Calculating the matrix A */

	if ((X.rows()!=Num_State_))
	{
		cout<<"The dimension of X in Calculate_A is wrong."<<endl;
		cout<<"Dimension of states is: "<<Num_State_<<endl;
		cout<<"Dimension of X "<<X.rows()<<endl;
		ERROR();
	}

	MatrixXd A; A.resize(Num_State_,Num_State_);A.setZero();

	if (Num_Com_>1)
	{
		Theta_=Calculate_Theta(X);
	}
	else
	{
		Theta_(Num_Com_-1)=1;
	}

	for (int i=0;i<Num_Com_;i++)
	{
		A=A+A_Matrix_[i]*Theta_(i);
	}


	return A;
}
VectorXd LPV::Calculate_Theta(VectorXd X)
{
	VectorXd Theta;Theta.resize(Num_Com_);Theta.setZero();

	for (int i=0;i<Num_Com_;i++)
	{
		Theta(i)=Prior_[i]*GaussianPDF(X,Mu_[i],Sigma_[i]);
	}
	double sum=Theta.sum();
	if (sum<1e-100)
	{
		for (int i=0;i<Num_Com_;i++)
		{
			Theta(i)=1.0/Num_Com_;
		}
	}
	else
	{
		Theta=Theta/sum;
	}

	return Theta;
}
MatrixXd LPV::readMatrix(const char *filename)
{
	int cols = 0, rows = 0;


	// Read numbers from file into buffer.
	ifstream infile;
	infile.open(filename);
	while (! infile.eof())
	{
		string line;
		getline(infile, line);

		int temp_cols = 0;
		stringstream stream(line);
		while(! stream.eof())
			stream >> buff[cols*rows+temp_cols++];

		if (temp_cols == 0)
			continue;

		if (cols == 0)
			cols = temp_cols;

		rows++;
	}

	infile.close();

	rows--;

	cout<<"row "<<rows<<" "<<cols<<endl;

	// Populate matrix with numbers.
	MatrixXd result(rows,cols);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result(i,j) = buff[ cols*i+j ];

	return result;
};

void LPV::ERROR()
{
	while(ros::ok())
	{

	}
}
double LPV::GaussianPDF(VectorXd x,VectorXd Mu,MatrixXd Sigma)
{

	double p;
	MatrixXd gfDiff;gfDiff.resize(1,Num_State_);
	MatrixXd gfDiff_T;gfDiff_T.resize(Num_State_,1);
	MatrixXd SigmaIIInv;SigmaIIInv.resize(Num_State_,Num_State_);
	double detSigmaII=0;
	MatrixXd gfDiffp;gfDiffp.resize(1,1);gfDiffp.setZero();

	detSigmaII=Sigma.determinant();
	SigmaIIInv=Sigma.inverse();
	if (detSigmaII<0)
	{
		detSigmaII=0;
	}
	gfDiff=(x - Mu).transpose();
	gfDiff_T=x - Mu;
	gfDiffp =gfDiff*SigmaIIInv* gfDiff_T;
	gfDiffp(0,0)=fabs(0.5*gfDiffp(0,0));
	p = exp(-gfDiffp(0,0)) / sqrt(pow(2.0*PI, Num_State_)*( detSigmaII +1e-50));
	return p;
}













void GMM::initialize(int Num_Com,int Num_state)
{

	/* Declare the number of the components
	 * and the dimension  of the state  */

	Num_Com_=Num_Com;
	Num_State_=Num_state;

	Prior_ =new double[Num_Com_];
	Mu_ =new VectorXd[Num_Com_];		for(int s=0; s<Num_Com_; s++ ){	Mu_[s].resize(Num_State_);	}
	Sigma_ =new MatrixXd[Num_Com_];		for(int s=0; s<Num_Com_; s++ ){	Sigma_[s].resize(Num_State_,Num_State_);	}


}
void GMM::initialize_GMM(const char  *path_prior_,const char  *path_mu_,const char  *path_sigma_, const char *path_threshold){

	/* Initialize scheduling parameters
	 *  path_prior_ is the path of the prior matrix
	 *	path_mu_ is the path of the mean matrix
	 *	path_sigma_ is the path of the covariance matrix */

	cout<<"Initializing the model of workspace"<<endl;

	MatrixXd fMatrix;
	if (is_file_exist(path_prior_))
	{
		fMatrix=readMatrix(path_prior_);
	}
	else
	{
		cout<<"The provided path does not exist."<<endl;
		cout<<"path_prior_GMM "<<endl;
		cout<<path_prior_<<endl;
		ERROR();
	}

	cout<<"fMatrix "<<endl;cout<<fMatrix<<endl;

	if ((fMatrix.cols()!=Num_Com_)||(fMatrix.rows()!=1))
	{
		cout<<"Initialization of Prior is wrong."<<endl;
		cout<<"Dimension of states is: "<<Num_Com_<<endl;
		cout<<"Dimension of states of Prior is: "<<fMatrix.cols()<<endl;
		ERROR();
	}

	cout<<"Prior"<<endl;
	for (int i=0; i<Num_Com_; i++)
	{
		Prior_[i]=fMatrix(0,i);
	}
	for (int i=0; i<Num_Com_; i++)
	{
		cout<<Prior_[i]<<endl;
	}

	fMatrix.setZero();
	if (is_file_exist(path_mu_))
	{
		fMatrix=readMatrix(path_mu_);
	}
	else
	{
		cout<<"The provided path does not exist."<<endl;
		cout<<"path_mu_GMM "<<endl;
		cout<<path_mu_<<endl;
		ERROR();
	}


	if ((fMatrix.cols()!=Num_Com_)||(fMatrix.rows()!=Num_State_))
	{
		cout<<"Initialization of Mean is wrong."<<endl;
		cout<<"Number of components is: "<<Num_Com_<<endl;
		cout<<"Dimension of states is: "<<Num_State_<<endl;
		cout<<"Dimension of states of Mean is: "<<fMatrix.rows()<<"*"<<fMatrix.cols()<<endl;
		ERROR();
	}

	for(int s=0; s<Num_Com_; s++ )
	{
		Mu_[s]=fMatrix.col(s);

	}

	for(int s=0; s<Num_Com_; s++ )
	{
		cout<<"Mu["<<s<<"]."<<endl;
		cout<<Mu_[s]<<endl;

	}

	fMatrix.resize(1,1);fMatrix.setZero();


	if (is_file_exist(path_sigma_))
	{
		fMatrix=readMatrix(path_sigma_);
	}
	else
	{
		cout<<"The provided path does not exist."<<endl;
		cout<<"path_sigma_GMM "<<endl;
		cout<<path_sigma_<<endl;
		ERROR();
	}


	if ((fMatrix.rows()!=Num_Com_*Num_State_)||(fMatrix.cols()!=Num_State_))
	{
		cout<<"Initialization of the covariance matrix is wrong."<<endl;
		cout<<"the covariance matrix : "<<endl;cout<<fMatrix<<endl;
		cout<<"Number of components is: "<<Num_Com_<<endl;
		cout<<"Dimension of states is: "<<Num_State_<<endl;
		cout<<"Dimension of states of the covariance matrix is: "<<fMatrix.rows()<<"*"<<fMatrix.cols()<<endl;
		ERROR();
	}

	int s=0;

	for(int i=0; i<Num_Com_; i++ ){
		for(int j=0; j<Num_State_; j++ ){
			Sigma_[i].row(j)=fMatrix.row(s);
			s++;
		}
	}

	for(int i=0; i<Num_Com_; i++ ){
		cout<<"Sigma["<<i<<"]."<<endl;
		cout<<Sigma_[i]<<endl;
	}
	fMatrix.resize(1,1);fMatrix.setZero();


	if (is_file_exist(path_threshold))
	{
		fMatrix=readMatrix(path_threshold);
	}
	else
	{
		cout<<"The provided path does not exist."<<endl;
		cout<<"path_threshold "<<endl;
		cout<<path_threshold<<endl;
		ERROR();
	}

	threshold=fMatrix(0,0);

	cout<<"The workspace threshold is "<<threshold<<endl;

}
double GMM::PDF(VectorXd X)
{
	VectorXd Theta;Theta.resize(Num_Com_);Theta.setZero();

	for (int i=0;i<Num_Com_;i++)
	{
		Theta(i)=Prior_[i]*GaussianPDF(X,Mu_[i],Sigma_[i]);
	}
	double sum=Theta.sum();


	return sum-threshold;
}
MatrixXd GMM::readMatrix(const char *filename)
{
	int cols = 0, rows = 0;

	// Read numbers from file into buffer.
	ifstream infile;
	infile.open(filename);

	while (! infile.eof())
	{
		string line;
		getline(infile, line);

		int temp_cols = 0;
		stringstream stream(line);
		while(! stream.eof())
			stream >> buff[cols*rows+temp_cols++];

		if (temp_cols == 0)
			continue;

		if (cols == 0)
			cols = temp_cols;

		rows++;
	}

	infile.close();

	rows--;

	// Populate matrix with numbers.
	MatrixXd result(rows,cols);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result(i,j) = buff[ cols*i+j ];

	return result;
};
void GMM::ERROR()
{
	while(ros::ok())
	{

	}
}
double GMM::GaussianPDF(VectorXd x,VectorXd Mu,MatrixXd Sigma)
{

	double p;
	MatrixXd gfDiff;gfDiff.resize(1,Num_State_);
	MatrixXd gfDiff_T;gfDiff_T.resize(Num_State_,1);
	MatrixXd SigmaIIInv;SigmaIIInv.resize(Num_State_,Num_State_);
	double detSigmaII=0;
	MatrixXd gfDiffp;gfDiffp.resize(1,1);gfDiffp.setZero();

	detSigmaII=Sigma.determinant();
	SigmaIIInv=Sigma.inverse();
	if (detSigmaII<0)
	{
		detSigmaII=0;
	}
	gfDiff=(x - Mu).transpose();
	gfDiff_T=x - Mu;
	gfDiffp =gfDiff*SigmaIIInv* gfDiff_T;
	gfDiffp(0,0)=fabs(0.5*gfDiffp(0,0));

	p = exp(-gfDiffp(0,0)) / sqrt(pow(2.0*PI, Num_State_)*( detSigmaII +1e-50));
	return p;
}



