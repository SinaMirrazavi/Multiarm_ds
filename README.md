# Multiarm_ds
[![Build Status](https://travis-ci.org/sinamr66/Multiarm_ds.svg?branch=master)](https://travis-ci.org/sinamr66/Multiarm_ds)

This repository includes the packages and instructions to generate motion for multi robots dual behavior system to coordinately reach a moving object or individually reach  independent targets. You can find the paper here: 

and the corresponding video here: 


# Dependences 

- Eigen http://eigen.tuxfamily.org/index.php?title=Main_Page

- motion_prediction https://github.com/epfl-lasa/Object-Trajectory-Prediction **(This package is dependent on Mathlib from form Robot-toolkit https://github.com/epfl-lasa/robot-toolkit)**


# Features:

- Multi arm motion generator.
  - Calculating the coordination parameter.
- Dual Behavior: Synchronous or Asynchronous.
  - Calculating the Synchronization parameter.
- Predicting the motion of the object and calculating the feasible intercept point.
  - It is possible to run this part of the class in the main program or run it independently.
 
# How to run 
### Initialization:
1- Initialize the Centralized dynamical system:
```
Initialize(Number of the robots, Number of the grabbing positions on the object, sample time, Number of end-effector state, A_V matix,Type of object motion prediction={Straight,Ballistic});
```
2- Initialize each robot Centralized dynamical system:
```
Initialize_robot(Index of the robot,Number of LPV Components, path to A matrix of the model of the robot's motion,  path to Priors of the model of the robot's motion, path to Mus, Path to Sigmas of the model of the robot's motion, Number of GMM Components of the model of the robot's workspace, Number of GMM states of the model of the robot's workspace, Path to priors of the model of the robot's workspace,Path to Mus to the model of the robot's workspace,Path Sigmas the model of the robot's workspace, Path to the feasibility threshold the model of the robot's workspace,Position of the base);
```
 **2-1 If you want to only use this code to predict the object motion and find the feasible intercept point, you only need to initialize the motion predication:** 
 ```
 Initialize_motion_predication_only(Number of the robots, Number of the grabbing positions on the object, sample time, Number of end-effector state,Type of object motion prediction={Straight,Ballistic});
 ```
3- Initialize the virtual object:
```
Initialize_the_virtual_object();
```
4- Set the initial position of the robot
```
  Set_the_initial_robot_state(Index of the robot,Initial position of the robot);
```
### In the loop:
5- Set the state of the robot
```
Set_the_robot_state(Index of the robot,State of the robot);
```
6- Set asynchronous target of the robot
```
Set_the_robot_first_primitive_desired_position(Index of the robot,State of the asynchronous target,D-State of the asynchronous target);
```
7- Set the state of the object
```
Set_the_object_state(State of the object,D-State of the object);
```
8- Set the state of the object for predication (I personally use filtered state for predication)
```
Set_the_object_state_for_prediction((State of the object,D-State of the object, Time);
```
Note: It is time! not the sample time.
9- Set state of the grabbing positions on the objects
```
Set_the_grabbing_state(Index of the robot,State of the grabbing position,State of the object );
```
10- Predict the object position 
```
predict_the_object_position();
```

11- Calculate the next desired motion for each robot
```
Update()
```
12- Get the desired state
```
Get_the_robot_state(Index of the robot, Desired state of the robot);
```
- **Note 1: Most of the variables are accessible by get_ functions see [multiarm_ds.h](https://github.com/sinamr66/Multiarm_ds/blob/master/include/multiarm_ds.h)**

- **Note 2: If you want to predict the object trajectory in a separate program, you need to Set_index_of_grabbing_posititon_(.)  and  Set_pos_of_grabbing_posititon_for_object_(.)  and Set_index_of_grabbing_posititon_(.) before Update(.) in the main program**

## Copyright
Please cite these papers if you are using this toolbox:
@article{mirrazavi2018unified,
  title={A unified framework for coordinated multi-arm motion planning},
  author={Mirrazavi Salehian, Seyed Sina and Figueroa, Nadia and Billard, Aude},
  journal={The International Journal of Robotics Research},
  pages={0278364918765952},
  publisher={SAGE Publications Sage UK: London, England}
}


For more information contact Sina Mirrazavi.
