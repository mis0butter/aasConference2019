% MAIN 
% Junette Hsin 
% Mohammad Ayoubi 

% clear; close all; 

%% To generate a new profile 

inputs                          % Creates initial, final, and sun vectors 
steering_profile                % Computes vel, accel, quat, torque profiles 
post_processing                 % Generates gyrostat unit sphere plot 

%% alpha0

load alpha0.mat 
post_processing                 % Generates gyrostat unit sphere plot 

%% alphaNot0 

load alphaNot0.mat 
post_processing                 % Generates gyrostat unit sphere plot 