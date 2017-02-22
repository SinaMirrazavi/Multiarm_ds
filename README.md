# Multiarm_ds
[![Build Status](https://travis-ci.org/sinamr66/Multiarm_ds.svg?branch=master)](https://travis-ci.org/sinamr66/Multiarm_ds)

This repository includes the packages and instructions to multi robots dual behavior coordinated motion to coordinately reach a moving object or individually reach a target. You can find the paper here: 

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
 
