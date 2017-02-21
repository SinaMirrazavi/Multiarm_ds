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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "eigen3/Eigen/Dense"
#include "ros/ros.h"


using namespace std;
using namespace Eigen;

#define MAXBUFSIZE  ((int) 1e5)




class LPV
{
public:
	void initialize(int Num_Com,int Num_state);
	void initialize_A(const char  *path_);
	void initialize_theta(const char  *path_prior_,const char  *path_mu_,const char  *path_sigma_);


	MatrixXd Calculate_A(VectorXd X);



private:

		void 	 	ERROR();
		double 		GaussianPDF(VectorXd x,VectorXd Mu,MatrixXd Sigma);
		VectorXd 	Calculate_Theta(VectorXd X);
		MatrixXd 	readMatrix(const char *filename);




		int 		Num_Com_;
		int 		Num_State_;
		VectorXd	Theta_;



		double 		*Prior_;
		VectorXd 	*Mu_;
		MatrixXd 	*Sigma_;

		MatrixXd 	*A_Matrix_;
		double 		buff[MAXBUFSIZE];



};




class GMM
{
public:
	void initialize(int Num_Com,int Num_state);
	void initialize_GMM(const char  *path_prior_,const char  *path_mu_,const char  *path_sigma_,const char *path_threshold);

	double 	PDF(VectorXd X);

private:

		void 	 	ERROR();
		double 		GaussianPDF(VectorXd x,VectorXd Mu,MatrixXd Sigma);

		MatrixXd 	readMatrix(const char *filename);




		int 		Num_Com_;
		int 		Num_State_;

		double 		*Prior_;
		VectorXd 	*Mu_;
		MatrixXd 	*Sigma_;

		double buff[MAXBUFSIZE];

		double threshold;




};


