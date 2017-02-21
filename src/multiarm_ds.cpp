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


#include "multiarm_ds.h"



void multiarm_ds::Initialize(int N_robots,int N_grabbing_pos, double dt, int N_state, MatrixXd A_V,ENUM_State_of_prediction Object_motion)
{
	/* Declare the number of the robots
	 * N_grabbing_pos is the number of the grabbing position.
	 * sample time (dt)
	 * the dimension  of the state
	 * and the gain matrix of the virtual object   */

	if (N_grabbing_pos>N_grabbing_pos)
	{
		cout<<"Number of the grabbing positions is more than the available robots."<<endl;
		ERROR();
	}
	Object_.order_of_grabbing_.resize(factorial(N_robots)/factorial(N_robots-N_grabbing_pos),N_grabbing_pos);
	Object_.order_of_grabbing_=PermGenerator(N_robots,N_robots);


	N_robots_=N_robots;
	N_state_=N_state;
	dt_=dt;
	Robots_=new S_Robot_ds[N_robots_];

	A_V_.resize(N_state_,N_state_);				A_V_.setZero();
	X_V_.resize(N_state_);						X_V_.setZero();
	DX_V_.resize(N_state_);						DX_V_.setZero();

	Object_.X_O_.resize(N_state_);				Object_.X_O_.setZero();
	Object_.X_O_First_.resize(N_state_);		Object_.X_O_First_.setZero();
	Object_.DX_O_.resize(N_state_);				Object_.DX_O_.setZero();
	Object_.X_I_C_.resize(N_state_);			Object_.X_I_C_.setZero();
	Object_.X_O_INTERCEPT_.resize(N_state_);	Object_.X_O_INTERCEPT_.setConstant(10);

	Vobject_.X_V_.resize(N_state);				Vobject_.X_V_.setZero();
	Vobject_.DX_V_INTERCEPT_.resize(N_state);	Vobject_.DX_V_INTERCEPT_.setZero();

	Object_.P_O_prediction_.resize((int)floor(Object_.MAX_PREDICTIONTIME_/(dt_))+1,3);
	for (int i=0;i<N_grabbing_pos;i++)
	{
		Vobject_.X_V_G_[i].resize(N_state);		Vobject_.X_V_G_[i].setZero();
		Vobject_.U_[i].resize(N_state);			Vobject_.U_[i].setZero();
		Vobject_.U_sum_.resize(N_state);		Vobject_.U_sum_.setZero();
		Object_.X_O_G_[i].resize(N_state);		Object_.X_O_G_[i].setZero();
		Object_.P_O_G_prediction_[i].resize((int)floor(Object_.MAX_PREDICTIONTIME_/(100*dt_))+1,3);
		Object_.P_O_G_prediction_[i].setZero();
	}
	for (int i=0;i<N_robots_;i++)
	{
		Robots_[i].X_.resize(N_state_);							Robots_[i].X_.setZero();
		Robots_[i].DX_.resize(N_state_);						Robots_[i].DX_.setZero();
		Robots_[i].X_d_.resize(N_state_);						Robots_[i].X_d_.setZero();
		Robots_[i].DX_d_.resize(N_state_);						Robots_[i].DX_d_.setZero();
		Robots_[i].X_F_P_.resize(N_state_);						Robots_[i].X_F_P_.setZero();
		Robots_[i].DX_F_P_.resize(N_state_);					Robots_[i].DX_F_P_.setZero();
		Robots_[i].ATX_.resize(N_state_);						Robots_[i].ATX_.setZero();
		Robots_[i].X_I_C_.resize(N_state_);						Robots_[i].X_I_C_.setZero();
		Robots_[i].Workspace_model_is_set_=false;
		Robots_[i].the_state_is_set_=false;
		Robots_[i].the_LPV_is_set_=false;
		Robots_[i].the_desired_state_is_set_=false;
	}

	if ((A_V.rows()!=N_state_)||(A_V.cols()!=N_state_))
	{
		cout<<"Initialization of "<<A_V<<" is wrong."<<endl;
		cout<<"A_V: "<<endl;cout<<A_V<<endl;
		ERROR();
	}
	A_V_=A_V;

	if (N_grabbing_pos>Max_Grabbing_state)
	{
		cout<<"The number of the grabbing state is higher than the Max. Just got to multiarm_ds.h and change Max_Grabbing_state"<<endl;
		ERROR();
	}
	Object_.N_grabbing_pos_=N_grabbing_pos;

	for (int i=0;i<Object_.N_grabbing_pos_;i++)
	{
		Object_.Grabbing_state_is_set_[i]=false;
		Vobject_.Grabbing_state_is_set[i]=false;
	}


	double gravity[3];gravity[0]=0;gravity[1]=0;
	Object_.Object_state_is_set_=false;
	Object_.Object_motion_=Object_motion;
	if (Object_.Object_motion_==Ballistic)
	{
		gravity[2]=9.81;
	}
	else
	{
		gravity[2]=0;
	}

	Object_.First_Object_state_is_set_=false;
	Object_.predict_ = new TrajectoryEstimation(gravity, dt_,30,10);
	cout<<"Object_.predict  is done"<<endl;
	Matrix3d handle;handle.setIdentity();
	Object_.predict_->setMomentOfInertia(handle);
	Object_.predict_->resetTrajectoryEstimation();
	The_catching_pos_is_found=false;
	Object_.Max_liklihood=0;
}
void multiarm_ds::Initialize_motion_predication_only(int N_robots,int N_grabbing_pos, double dt, int N_state,ENUM_State_of_prediction Object_motion)
{
	/* If you want to only predict the motion of the object, use this!
	 * Declare the number of the robots
	 * N_grabbing_pos is the number of the grabbing position.
	 * sample time (dt)
	 * the dimension  of the state  */
	if (N_grabbing_pos>N_grabbing_pos)
	{
		cout<<"Number of the grabbing positions is more than the available robots."<<endl;
		ERROR();
	}
	Object_.order_of_grabbing_.resize(factorial(N_robots)/factorial(N_robots-N_grabbing_pos),N_grabbing_pos);
	Object_.order_of_grabbing_=PermGenerator(N_robots,N_robots);


	N_robots_=N_robots;
	N_state_=N_state;
	dt_=dt;
	Robots_=new S_Robot_ds[N_robots_];

	A_V_.resize(N_state_,N_state_);				A_V_.setZero();
	X_V_.resize(N_state_);						X_V_.setZero();
	DX_V_.resize(N_state_);						DX_V_.setZero();

	Object_.X_O_.resize(N_state_);				Object_.X_O_.setZero();
	Object_.X_O_First_.resize(N_state_);		Object_.X_O_First_.setZero();
	Object_.DX_O_.resize(N_state_);				Object_.DX_O_.setZero();
	Object_.X_I_C_.resize(N_state_);			Object_.X_I_C_.setZero();
	Object_.X_O_INTERCEPT_.resize(N_state_);	Object_.X_O_INTERCEPT_.setConstant(10);

	Vobject_.X_V_.resize(N_state);				Vobject_.X_V_.setZero();
	Vobject_.DX_V_INTERCEPT_.resize(N_state);	Vobject_.DX_V_INTERCEPT_.setZero();

	Object_.P_O_prediction_.resize((int)floor(Object_.MAX_PREDICTIONTIME_/(dt_))+1,3);
	for (int i=0;i<N_grabbing_pos;i++)
	{
		Vobject_.X_V_G_[i].resize(N_state);		Vobject_.X_V_G_[i].setZero();
		Vobject_.U_[i].resize(N_state);			Vobject_.U_[i].setZero();
		Vobject_.U_sum_.resize(N_state);		Vobject_.U_sum_.setZero();
		Object_.X_O_G_[i].resize(N_state);		Object_.X_O_G_[i].setZero();
		Object_.P_O_G_prediction_[i].resize((int)floor(Object_.MAX_PREDICTIONTIME_/(dt_))+1,3);
		Object_.P_O_G_prediction_[i].setZero();
	}
	for (int i=0;i<N_robots_;i++)
	{
		Robots_[i].X_.resize(N_state_);							Robots_[i].X_.setZero();
		Robots_[i].DX_.resize(N_state_);						Robots_[i].DX_.setZero();
		Robots_[i].X_d_.resize(N_state_);						Robots_[i].X_d_.setZero();
		Robots_[i].DX_d_.resize(N_state_);						Robots_[i].DX_d_.setZero();
		Robots_[i].X_F_P_.resize(N_state_);						Robots_[i].X_F_P_.setZero();
		Robots_[i].DX_F_P_.resize(N_state_);					Robots_[i].DX_F_P_.setZero();
		Robots_[i].ATX_.resize(N_state_);						Robots_[i].ATX_.setZero();
		Robots_[i].X_I_C_.resize(N_state_);						Robots_[i].X_I_C_.setZero();
		Robots_[i].Workspace_model_is_set_=false;
		Robots_[i].the_state_is_set_=false;
		Robots_[i].the_LPV_is_set_=false;
		Robots_[i].the_desired_state_is_set_=false;
	}

	if (N_grabbing_pos>Max_Grabbing_state)
	{
		cout<<"The number of the grabbing state is higher than the Max. Just got to multiarm_ds.h and change Max_Grabbing_state"<<endl;
		ERROR();
	}
	Object_.N_grabbing_pos_=N_grabbing_pos;

	for (int i=0;i<Object_.N_grabbing_pos_;i++)
	{
		Object_.Grabbing_state_is_set_[i]=false;
		Vobject_.Grabbing_state_is_set[i]=false;
	}


	double gravity[3];gravity[0]=0;gravity[1]=0;
	Object_.Object_state_is_set_=false;
	Object_.Object_motion_=Object_motion;
	if (Object_.Object_motion_==Ballistic)
	{
		gravity[2]=9.81;
	}
	else
	{
		gravity[2]=0;
	}

	Object_.First_Object_state_is_set_=false;
	Object_.predict_ = new TrajectoryEstimation(gravity, dt_,100,10);
	cout<<"Object_.predict  is done"<<endl;
	Matrix3d handle;handle.setIdentity();
	Object_.predict_->setMomentOfInertia(handle);
	Object_.predict_->resetTrajectoryEstimation();
	The_catching_pos_is_found=false;
	Object_.Max_liklihood=0;

}

void multiarm_ds::Initialize_robot(int index,int Num_LPV_Com, const char  *path_A_LPV, const char  *path_prior_LPV,const char  *path_mu_LPV,const char  *path_sigma_LPV,int Num_GMM_Com, int Num_GMM_state, const char  *path_prior_GMM,const char  *path_mu_GMM,const char  *path_sigma_GMM, const char *path_threshold,Vector3d X_Base)
{
	/*Declare the parameters of   index th robot.
	 * Num_LPV_Com Number of components of LPV
	 * path_A_LPV is the path to A matrices
	 * path_prior_LPV is the prior of Theta function
	 * path_mu_LPV is the Mean of Theta function
	 * path_sigma_LPV is the Sigma of Theta function
	 * Num_GMM_Com	 is the number of the components of the robot's workspace model
	 * Num_GMM_state is the number of the states of the robot's workspace model
	 * Num_GMM_state is the number of the states of the robot's workspace model
	 * path_prior_GMM is the number of the states of the robot's workspace model
	 * path_mu_GMM is the number of the states of the robot's workspace model
	 * path_sigma_GMM is the number of the states of the robot's workspace model
	 * X_Base is the position of the base of the robot*/



	if (index>N_robots_-1)
	{
		cout<<"Initialization of "<<index<<"th robot is wrong."<<endl;
		cout<<"index "<<index<<" Max robot Number "<<N_robots_-1<<endl;
		ERROR();
	}


	if (Robots_[index].Workspace_model_is_set_||Robots_[index].the_LPV_is_set_)
	{
		cout<<"Initialization of "<<index<<"th robot is wrong."<<endl;
		cout<<"This robot is already being initialized"<<endl;
		ERROR();
	}
	if (X_Base.rows()!=3)
	{
		cout<<"Initialization of "<<index<<"th robot is wrong."<<endl;
		cout<<"The base of the robot is wrong. The size of X_Base is "<<X_Base.rows()<<endl;
		ERROR();
	}



	Robots_[index].Workspace_model_.initialize(Num_GMM_Com,Num_GMM_state);
	Robots_[index].Dynamic_.initialize(Num_LPV_Com,N_state_);

	Robots_[index].Workspace_model_.initialize_GMM(path_prior_GMM,path_mu_GMM,path_sigma_GMM,path_threshold);
	Robots_[index].Dynamic_.initialize_theta(path_prior_LPV,path_mu_LPV,path_sigma_LPV);
	Robots_[index].Dynamic_.initialize_A(path_A_LPV);
	Robots_[index].X_Base_=X_Base;

	cout<<"The base of "<<index<<"th robot is "<<endl;
	cout<<Robots_[index].X_Base_<<endl;


	Robots_[index].Workspace_model_is_set_=true;
	Robots_[index].the_LPV_is_set_=true;
	Robots_[index].tau_=0.0001;
	Robots_[index].Dtau_=0;
}
void multiarm_ds::Set_the_robot_state(int index,VectorXd X)
{
	/* Setting the current state of the  index th robot
	 * X is the state of the end-effector with respect to the world-frame
	 * 								*/
	if ((Robots_[index].X_.rows()!=X.rows()))
	{
		cout<<"The state dimension of "<<index<<"th robot is wrong."<<endl;
		cout<<"The dimension of X is "<<X.rows()<<endl;
		cout<<"The dimension of robot is "<<Robots_[index].X_.rows()<<endl;
		ERROR();
	}
	if (Robots_[index].the_state_is_set_==true)
	{
		cout<<"States of "<<index<<"th robot is already being set."<<endl;
		ERROR();
	}

	Robots_[index].X_=X;
	Robots_[index].the_state_is_set_=true;
}
void multiarm_ds::Set_the_initial_robot_state(int index,Vector3d X)
{
	/* Setting the initial state of the  index th robot
	 * X is the state of the end-effector with respect to the world-frame
	 * 								*/

	Robots_[index].X_Initial_pose_=X;
	cout<<"The  initial position of "<<index<<"th robot is:"<<endl<<Robots_[index].X_Initial_pose_<<endl;
}


void multiarm_ds::Set_the_robot_first_primitive_desired_position(int index,VectorXd X,VectorXd DX)
{

	/* Setting the desired state of the first primitive of the  index th robot
	 * X is the desired state of the end-effector with respect to the world-frame
	 * DX is the D-desired state of the end-effector
	 * 	DOES NOT NEED TO CALL AT EACH UPDATE LOOP		*/

	if ((Robots_[index].X_F_P_.rows()!=X.rows())||(Robots_[index].DX_F_P_.rows()!=DX.rows()))
	{
		cout<<"The dimension of the desired state of "<<index<<"th robot is wrong."<<endl;
		cout<<"The dimension of X is "<<X.rows()<<" and "<<DX.rows()<<endl;
		cout<<"The dimension of robot is "<<Robots_[index].X_F_P_.rows()<<endl;
		ERROR();
	}


	Robots_[index].X_F_P_=X;
	Robots_[index].DX_F_P_=DX;
	Robots_[index].the_desired_state_is_set_=true;
}
void multiarm_ds::Set_the_object_state(VectorXd X,VectorXd DX)
{
	/* Setting the current state of the object
	 * X is the state of the object with respect to the world-frame
	 * DX is the  derivative of the state of the object		*/

	if (Object_.Object_state_is_set_==true)
	{
		cout<<"States of the object is already being set."<<endl;
		ERROR();
	}

	if ((X.rows()!=N_state_)&&(DX.rows()!=N_state_))
	{
		cout<<"The dimension of the state of the object is wrong"<<endl;
		ERROR();
	}

	Object_.X_O_=X;
	Object_.DX_O_=DX;

	Object_.Object_state_is_set_=true;
}
void multiarm_ds::Set_the_object_state_for_prediction(VectorXd X,VectorXd X_filtered, double time)
{
	/* Setting the current state of the object
	 * X is the state of the object with respect to the world-frame
	 * DX is the  derivative of the state of the object		*/

	if ((X_filtered.rows()!=3))
	{
		cout<<"The dimension of the state of the object is wrong"<<endl;
		ERROR();
	}


	Object_.X_O_.block(0,0,3,1)=X;
	if (Object_.First_Object_state_is_set_==false)
	{
		Object_.First_Object_state_is_set_=true;
		Object_.X_O_First_=Object_.X_O_.block(0,0,3,1);
		cout<<"The first position of the object is saved:"<<endl;
		cout<<Object_.X_O_First_<<endl;
	}
	else if ((Object_.First_Object_state_is_set_==true)&&(fabs(Object_.X_O_(0)-Object_.X_O_First_(0))>0.4))
	{
		Object_.predict_->setNewObservationset(time, X_filtered);
	}
	Object_.Object_state_is_set_=true;
}


void multiarm_ds::Set_index_of_grabbing_posititon_(int index_of_robot,int indext_of_grabbing_pos)
{
	Robots_[index_of_robot].index_of_grabbing_posititon_=indext_of_grabbing_pos;
}

void multiarm_ds::Set_the_grabbing_state(int index,VectorXd X_G)
{
	/* Setting the grabbing states on the object with respect to the universal frame (the same frame as the object and etc...
	 * X is the state of index th grabbing states with respect to the world-frame
	 * Notice: you just need to set it once! and you need to set it exactly after the object state is set!		*/

	if (Object_.Object_state_is_set_==false)
	{
		cout<<"States of the object has not been set. yet"<<endl;
		ERROR();
	}

	if (Object_.Grabbing_state_is_set_[index]==true)
	{
		cout<<"States of the "<<index<<" th grabbing position is already being set."<<endl;
		ERROR();
	}

	Set_the_grabbing_state(index,X_G,Object_.X_O_);

}
void multiarm_ds::Set_the_grabbing_state(int index,VectorXd X_G,VectorXd X_O)
{
	/* Setting the grabbing states on the object with respect to the universal frame (the same frame as the object and etc...
	 * X is the state of index th grabbing states with respect to the world-frame
	 * Notice: you just need to set it once! and you need to set it exactly after the object state is set!		*/

	if (Object_.Object_state_is_set_==false)
	{
		cout<<"States of the object has not been set. yet"<<endl;
		ERROR();
	}

	if (Object_.Grabbing_state_is_set_[index]==true)
	{
		cout<<"States of the "<<index<<" th grabbing position is already being set."<<endl;
		ERROR();
	}

	Object_.X_O_G_[index]=X_G-X_O;
	Object_.Grabbing_state_is_set_[index]=true;

}
void multiarm_ds::Set_pos_of_grabbing_posititon_for_object_(bool catching_pos_is_found,double likelihood, Vector3d X_I_C)
{
	The_catching_pos_is_found=catching_pos_is_found;
	Object_.Max_liklihood=likelihood;
	Object_.X_I_C_.block(0,0,3,1)=X_I_C;

	//	cout<<"Object_.Max_liklihood "<<Object_.Max_liklihood<<endl;
	//	cout<<"Object_.X_I_C_.block(0,0,3,1)"<<Object_.X_I_C_.block(0,0,3,1)<<endl;
}
void multiarm_ds::Set_index_of_grabbing_posititon_(int index_of_robot, int index_of_grabbing, Vector3d X_I_C)
{

	/* Setting the desired grabbing position of the object for index_of_robot th robot with respect to the world frame
	 * X_I_C is the desired grabbing position in world frame
	 * index_of_grabbing the index of the desired grabbing position	*/

	Robots_[index_of_robot].index_of_grabbing_posititon_=index_of_grabbing;
	Robots_[index_of_robot].X_I_C_.block(0,0,3,1)=X_I_C;


	//	cout<<"Robot "<<index_of_robot<<" "<<Robots_[index_of_robot].index_of_grabbing_posititon_<<endl;
	//	cout<<"Robots_[index_of_robot].X_I_C_ "<<Robots_[index_of_robot].X_I_C_.block(0,0,3,1)<<endl;

}
void multiarm_ds::Initialize_the_virtual_object()
{

	/* Initializing the state of the virtual object
	 * This one needs to be called once, when the states of all the grabbing positions and the robots are set. 	*/
	if (Object_.Object_state_is_set_==false)
	{
		cout<<"States of the object has not been set. yet"<<endl;
		ERROR();
	}

	for (int i=0;i<Object_.N_grabbing_pos_;i++)
	{
		if (Object_.Grabbing_state_is_set_[i]==false)
		{
			cout<<"States of the "<<i<<" th grabbing position is not being set yet."<<endl;
			ERROR();
		}
	}

	for (int i=0;i<N_robots_;i++)
	{
		if (Robots_[i].the_state_is_set_==false)
		{
			cout<<"States of the "<<i<<" th robot is not being set yet."<<endl;
			ERROR();
		}
	}


	Vobject_.N_grabbing_pos=Object_.N_grabbing_pos_;


	for (int i=0;i<Object_.N_grabbing_pos_;i++)
	{
		Vobject_.X_V_G_[i]=Object_.X_O_G_[i];
		Vobject_.Grabbing_state_is_set[i]=true;
		cout<<"The distance between the object and the "<<i<<"th grabbing position is"<<endl;cout<<Vobject_.X_V_G_[i]<<endl;
	}

	VectorXd Handle; Handle.resize(N_state_); Handle.setZero();

	for (int i=0;i<N_robots_;i++)
	{
		Handle=Handle+Robots_[i].X_;
	}
	Handle=Handle/N_robots_;

	Vobject_.X_V_=Handle;
	cout<<"The position of the virtual object is"<<endl;cout<<Vobject_.X_V_<<endl;
}
void multiarm_ds::calculate_coordination_allocation()
{

	/* Calculate the coordination allocation parameters.
	 * tau
	 * Dtau 	*/


	Vobject_.tau_sum_=0;
	for (int i=0;i<N_robots_;i++)
	{
		MatrixXd handle(1,1);

		Robots_[i].M[0]=Object_.DX_O_(0);
		Robots_[i].M[1]=fabs(sgn(Robots_[i].index_of_grabbing_posititon_+1))*fabs(Object_.Max_liklihood)*min(sgn(Object_.Max_liklihood),sgn(Robots_[i].index_of_grabbing_posititon_+1));
		Robots_[i].Dtau_=10*Robots_[i].tau_*(1-Robots_[i].tau_)*Robots_[i].M[0]/((Object_.X_O_.block(0,0,3,1)-Vobject_.X_V_.block(0,0,3,1)).norm());
		Robots_[i].tau_=Robots_[i].tau_+Robots_[i].Dtau_*dt_;
		if (Robots_[i].tau_<=0.001)
		{
			Robots_[i].tau_=0.001;
			Robots_[i].Dtau_=0.0;
		}
		else if (Robots_[i].tau_>0.99)
		{

			Robots_[i].tau_=0.99;
			Robots_[i].Dtau_=0.0;
		}
		Vobject_.tau_sum_=Vobject_.tau_sum_+Robots_[i].tau_;
	}

}
void multiarm_ds::calculate_coordination_parameter()
{
	/* Calculate Gamma and Dgamma*/
	//	Vobject_.Dgamma_=0.01/(Object_.X_O_INTERCEPT_.block(0,0,3,1).norm()+1e-32);
	MatrixXd handle_exp=Object_.X_O_INTERCEPT_.block(0,0,3,1).transpose()*Object_.X_O_INTERCEPT_.block(0,0,3,1);
	MatrixXd Handle_M=(Object_.DX_O_.block(0,0,3,1).transpose()*Object_.X_O_INTERCEPT_.block(0,0,3,1));
	Vobject_.Dgamma_=-2*Handle_M(0,0)/exp(handle_exp(0,0));
	Vobject_.gamma_=1/exp(0.1*handle_exp(0,0));
	//Vobject_.gamma_=Vobject_.gamma_+Vobject_.Dgamma_*dt_;
	if (Vobject_.gamma_>=1)
	{
		Vobject_.gamma_=1;
		Vobject_.Dgamma_=0;
	}
	if (Vobject_.gamma_<0)
	{
		Vobject_.gamma_=0;
		Vobject_.Dgamma_=0;
	}
	//	cout<<"Vobject_.gamma "<<Vobject_.gamma_<<" "<<Vobject_.Dgamma_<<endl;

}
void multiarm_ds::calculate_robot_next_state()
{
	for (int i=0;i<N_robots_;i++)
	{
		Robots_[i].DX_=Robots_[i].tau_*Robots_[i].DX_d_+Robots_[i].ATX_;
	}
}
void multiarm_ds::calculate_ATX()
{
	for (int i=0;i<N_robots_;i++)
	{
		Robots_[i].ATX_=Robots_[i].Dynamic_.Calculate_A(Robots_[i].X_)*(Robots_[i].X_-Robots_[i].X_F_P_-Robots_[i].tau_*(Robots_[i].X_d_-Robots_[i].X_F_P_))+Robots_[i].Dtau_*(Robots_[i].X_d_-Robots_[i].X_F_P_);
	}
}
void multiarm_ds::calculate_u()
{
	Vobject_.U_sum_.setZero();

	for (int i=0;i<N_robots_;i++)
	{
		//		Vobject_.U[i]=Robots_[i].DX-Robots_[i].ATX;
		Vobject_.U_[i]=-Robots_[i].ATX_;
		Vobject_.U_sum_=Vobject_.U_sum_+Robots_[i].tau_*Vobject_.U_[i];
	}
}


void multiarm_ds::CalculateCatchingprobability(int nFrame,int i_robot,int i_object)
{
	MatrixXd r_P_o;r_P_o.resize(nFrame,3);
	/*#pragma omp parallel num_threads(8)
	{
#pragma omp for*/
	for (int i=0;i<nFrame;i++)
	{
		r_P_o.row(i)=Object_.P_O_G_prediction_[i_object].row(i)-Robots_[i_robot].X_Base_.transpose();
		if( ( (r_P_o(i,0)>minPos[0] ) &&  (r_P_o(i,0)<maxPos[0] ) ) &&
				( (r_P_o(i,1)>minPos[1] ) &&  (r_P_o(i,1)<maxPos[1] ) ) &&
				( (r_P_o(i,2)>minPos[2] ) &&  (r_P_o(i,2)<maxPos[2] ) ))
		{
			Robots_[i_robot].Probability_of_catching_(i,i_object)=Robots_[i_robot].Workspace_model_.PDF(r_P_o.row(i).transpose());
			//			cout<<"nFrame "<<nFrame<<" "<<i<<" "<<Robots_[i_robot].Probability_of_catching_(i,i_object)<<endl;
			if (Robots_[i_robot].Probability_of_catching_(i,i_object)<0)
			{
				Robots_[i_robot].Probability_of_catching_(i,i_object)=0;
			}
		}

	}
	//}

}

void multiarm_ds::predict_the_object_position()
{
	if ((Object_.predict_->mReadyToPredict)&&(Object_.X_O_(0)<4*minPos[0]))
	{
		int nFrame =0;
		nFrame = Object_.predict_->PredictNextPosVel(Object_.MAX_PREDICTIONTIME_,dt_, Object_.P_O_prediction_,true,0,maxPos[0]);
		//	cout<<"frame "<<nFrame<<" a "<<Object_.P_O_prediction_.cols()<<" "<<Object_.P_O_prediction_.rows()<<endl;
		//	cout<<"frame "<<nFrame<<endl;

		if 	(nFrame==-1)
		{
			for (int i=0;i<N_robots_;i++)
			{
				Object_.Max_liklihood=0.0;
			}
		}
		else
		{
			for (int i=0;i<Object_.N_grabbing_pos_;i++)
			{
				for (int j=0;j<nFrame;j++)
				{
					Object_.P_O_G_prediction_[i].row(j)=Object_.P_O_prediction_.row(j)+Object_.X_O_G_[i].block(0,0,3,1).transpose();
				}
				//	cout<<"Object_.P_O_G_prediction_[i] "<<i<<" "<<Object_.P_O_G_prediction_[i]<<endl;
			}
			//}

			for (int i=0;i<N_robots_;i++)
			{
				Robots_[i].Probability_of_catching_.resize(nFrame,Object_.N_grabbing_pos_);Robots_[i].Probability_of_catching_.setZero();
				for (int j=0;j<Object_.N_grabbing_pos_;j++)
				{
					CalculateCatchingprobability(nFrame,i,j);
				}
				//				cout<<"Robots_[i].Probability_of_catching "<<i<<endl<<Robots_[i].Probability_of_catching_<<endl;
			}

			Object_.prob_order_of_grabbing_.resize(nFrame,Object_.order_of_grabbing_.rows());Object_.prob_order_of_grabbing_.setOnes();
			int pointer_robot=-1;
			for (int i=0;i<Object_.order_of_grabbing_.rows();i++)
			{
				for (int j=0;j<Object_.order_of_grabbing_.cols();j++)
				{
					pointer_robot=Object_.order_of_grabbing_(i,j);
					Object_.prob_order_of_grabbing_.col(i)=Object_.prob_order_of_grabbing_.col(i).cwiseProduct(Robots_[pointer_robot].Probability_of_catching_.col(j));
				}
			}

			//	cout<<"Object_.prob_order_of_grabbing "<<endl<<Object_.prob_order_of_grabbing_<<endl;
			int index_row=0;
			int index_column=0;
			int handle=0;
			double handle_distance=10000;
			double handle_distance1=0;
			Object_.Max_liklihood=Object_.prob_order_of_grabbing_.maxCoeff(&index_row ,& index_column);
			//	cout<<"max_likihood "<<Object_.Max_liklihood<<" index_row "<<index_row<<" index_column "<<index_column<<endl;
			handle=0;

			for (int i=0;i<N_robots_;i++)
			{
				Robots_[i].index_of_grabbing_posititon_=-2;
			}
			if (Object_.Max_liklihood>0.00001)
			{
				for (int i=0;i<Object_.N_grabbing_pos_;i++)
				{
					handle=Object_.order_of_grabbing_(index_column,i);
					Robots_[handle].index_of_grabbing_posititon_=i;
				}
				for (int i=0;i<Object_.prob_order_of_grabbing_.rows();i++)
				{
					if (Object_.prob_order_of_grabbing_(i,index_column)>0)
					{
						handle_distance1=0;
						if (The_catching_pos_is_found)
						{
							for (int j=0;j<N_robots_;j++)
							{
								handle_distance1=handle_distance1+(Robots_[j].X_I_C_.block(0,0,3,1)-Object_.P_O_G_prediction_[Robots_[j].index_of_grabbing_posititon_].row(i).transpose()).norm();
							}
						}
						else
						{
							for (int j=0;j<N_robots_;j++)
							{
								handle_distance1=handle_distance1+(Robots_[j].X_Initial_pose_-Object_.P_O_G_prediction_[Robots_[j].index_of_grabbing_posititon_].row(i).transpose()).norm();
							}
						}
						if (handle_distance>handle_distance1)
						{
							index_row=i;
							handle_distance=handle_distance1;
						}

					}
				}
				cout<<"Catching position likelihood "<<Object_.prob_order_of_grabbing_(index_row,index_column)<<" index_row "<<index_row<<" index_column "<<index_column<<endl;
				Object_.X_I_C_.block(0,0,3,1)=Object_.P_O_prediction_.row(index_row).transpose();
				Object_.index_row=index_row;
				Object_.index_column=index_column;
				for (int i=0;i<Object_.N_grabbing_pos_;i++)
				{
					handle=Object_.order_of_grabbing_(index_column,i);
					Robots_[handle].index_of_grabbing_posititon_=i;
					Robots_[handle].X_I_C_.block(0,0,3,1)=Object_.P_O_G_prediction_[i].row(index_row).transpose();
				}
				for (int i=0;i<N_robots_;i++)
				{
					cout<<"The "<<i<<" th robot is assigned to grab "<<Robots_[i].index_of_grabbing_posititon_<<" the catching position. ";
				}
				cout<<""<<endl;
				The_catching_pos_is_found=true;
			}
			else
			{
				Object_.Max_liklihood=-1;
			}
		}

	}
	else if ( ((Object_.X_O_(0)>2*minPos[0]) &&  (Object_.X_O_(0)<maxPos[0] ) ) &&
			( (Object_.X_O_(2)>minPos[2] ) &&  (Object_.X_O_(2)<maxPos[2] ) ))
	{
		Object_.Max_liklihood=10;
		int handle=0;
		for (int i=0;i<Object_.N_grabbing_pos_;i++)
		{
			handle=Object_.order_of_grabbing_(Object_.index_column,i);
			Robots_[handle].index_of_grabbing_posititon_=i;
			Robots_[handle].X_I_C_.block(0,0,3,1)=Object_.X_O_.block(0,0,3,1)+Object_.X_O_G_[i].block(0,0,3,1);
		}
		Object_.X_I_C_.block(0,0,3,1)=Object_.X_O_.block(0,0,3,1);

	}

}


void multiarm_ds::assign_the_robots()
{
	for (int i=0;i<N_robots_;i++)
	{
		if (Robots_[i].index_of_grabbing_posititon_>-1)
		{
			Robots_[i].DX_d_=Vobject_.DX_V_INTERCEPT_;
			Robots_[i].X_d_=Vobject_.X_V_G_[Robots_[i].index_of_grabbing_posititon_]+Vobject_.X_V_;
		}
	}
}



void multiarm_ds::Update()
{


	if (!everythingisreceived())
	{
		cout<<"You forgot to set sth! and you called Update"<<endl;
		ERROR();
	}
	if (The_catching_pos_is_found)
	{
		calculate_coordination_parameter();
		calculate_coordination_allocation();
		assign_the_robots();
		calculate_ATX();
		calculate_u();


		Object_.X_O_INTERCEPT_=Object_.X_O_-Object_.X_I_C_;
		Vobject_.X_V_INTERCEPT_=Vobject_.X_V_-Object_.X_I_C_;

		Vobject_.DX_V_INTERCEPT_=(Vobject_.gamma_*Object_.DX_O_+Vobject_.Dgamma_*Object_.X_O_INTERCEPT_
				+A_V_*(Vobject_.X_V_INTERCEPT_-Vobject_.gamma_*Object_.X_O_INTERCEPT_)
				+Vobject_.U_sum_)/(1+Vobject_.tau_sum_);

		/*		Vobject_.DX_V_INTERCEPT_=(Vobject_.Dgamma_*Object_.X_O_INTERCEPT_
					+A_V_*(Vobject_.X_V_INTERCEPT_-Vobject_.gamma_*Object_.X_O_INTERCEPT_)
					+Vobject_.U_sum_)/(1+Vobject_.tau_sum_);*/


		calculate_robot_next_state();



		Vobject_.X_V_=Vobject_.X_V_+Vobject_.DX_V_INTERCEPT_*dt_;

		for (int i=0;i<N_robots_;i++)
		{
			Robots_[i].X_=Robots_[i].X_+Robots_[i].DX_*dt_;
		}

	}


	restart_everything();
}


bool multiarm_ds::Get_prediction_state()
{
	return Object_.predict_->mReadyToPredict;
}

bool multiarm_ds::Get_catching_state()
{
	return The_catching_pos_is_found;
}
void multiarm_ds::Get_the_coordination_allocation(int index, double& x)
{
	x=Robots_[index].tau_;
}
void multiarm_ds::Get_the_coordination_parameter(double& x)
{
	x=Vobject_.gamma_;
}
void multiarm_ds::Get_predict_the_object_position(int index, MatrixXd& X)
{
	//X.resize(Object_.P_O_G_prediction_[index].rows(),Object_.P_O_G_prediction_[index].cols());
	X=Object_.P_O_G_prediction_[index];
}
void multiarm_ds::Get_the_robot_state(int index, VectorXd& X)
{

	/* Getting the desired state of index th robot
	 *
	 * X is the desired state of the robot	*/

	X=Robots_[index].X_;
}
void multiarm_ds::Get_Virtual_state(VectorXd & X)
{
	/* Getting the current state of the virtual object
	 * X is the state of the virtual object	*/

	X=Vobject_.X_V_;
}
void multiarm_ds::Get_the_grabbing_state(int index, VectorXd & X)
{
	/* Getting the current state of index th grabbing position on the virtual object with respect to the world frame
	 * X is the state of the virtual object	*/

	X=Vobject_.X_V_G_[index]+Vobject_.X_V_;
}
bool multiarm_ds::Get_pos_of_grabbing_posititon_for_object_(double& likelihood, Vector3d& X_I_C)
{
	/* Getting the desired grabbing position of the virtual object
	 * X_I_C is the desired grabbing position in world frame
	 * likelihood is the related to change to grab!	*/
	if (The_catching_pos_is_found)
	{
		X_I_C=Object_.X_I_C_.block(0,0,3,1);
		likelihood=Object_.Max_liklihood;
	}
	return The_catching_pos_is_found;
}
void multiarm_ds::Get_index_of_grabbing_posititon_(int index_of_robot, int& index_of_grabbing, Vector3d& X_I_C)
{

	/* Getting the desired grabbing position of the object for index_of_robot th robot with respect to the world frame
	 * X_I_C is the desired grabbing position in world frame
	 * index_of_grabbing the index of the desired grabbing position	*/


	if (The_catching_pos_is_found)
	{
		index_of_grabbing=Robots_[index_of_robot].index_of_grabbing_posititon_;
		X_I_C=Robots_[index_of_robot].X_I_C_.block(0,0,3,1);
	}


}


/*void multiarm_ds::Get_the_desired_intercept_state(int index, VectorXd & X)
{
	 Getting the desired intercept point of index th robot
 * X is the desired intercept point

	X=Robots_[index].X_I_C_.block(0,0,3,1);
}*/




















void multiarm_ds::ERROR()
{
	while(ros::ok())
	{

	}
}
void multiarm_ds::restart_everything()
{
	for(int i=0;i<N_robots_;i++)
	{
		Robots_[i].the_state_is_set_=false;
	}
	Object_.Object_state_is_set_=false;

}
bool multiarm_ds::everythingisreceived()
{
	bool flag=true;

	for (int i=0;i<N_robots_;i++)
	{
		if ((Robots_[i].the_state_is_set_==false)||(Robots_[i].Workspace_model_is_set_==false)||(Robots_[i].the_LPV_is_set_==false)||(Robots_[i].the_desired_state_is_set_==false))
		{
			flag=false;
			cout<<"Sth from robot is missing"<<endl;
			cout<<"Robots_[i].the_state_is_set "<<Robots_[i].the_state_is_set_<<endl;
			cout<<"Robots_[i].Workspace_model_is_set_ "<<Robots_[i].Workspace_model_is_set_<<endl;
			cout<<"Robots_[i].the_LPV_is_set_ "<<Robots_[i].the_LPV_is_set_<<endl;
			cout<<"Robots_[i].the_desired_state_is_set_ "<<Robots_[i].the_desired_state_is_set_<<endl;
		}
	}
	if (Object_.Object_state_is_set_==false)
	{
		cout<<"Object_state_is_set_ is missing"<<endl;
		flag=false;
	}
	for (int i=0;i<Object_.N_grabbing_pos_;i++)
	{
		if ((Object_.Grabbing_state_is_set_[i]==false)||(Vobject_.Grabbing_state_is_set[i]==false))
		{
			cout<<"Sth from grabbing state is missing"<<endl;
			cout<<"Object_.Grabbing_state_is_set_ "<<Object_.Grabbing_state_is_set_[i]<<endl;
			cout<<"Vobject_.Grabbing_state_is_set "<<Vobject_.Grabbing_state_is_set[i]<<endl;
			flag=false;
		}
	}

	return flag;
}
