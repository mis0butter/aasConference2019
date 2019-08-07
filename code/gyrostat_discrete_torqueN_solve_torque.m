function [time_out, q_out, w_out, torque_out] = gyrostat_discrete_torqueN_solve_torque(dt, t_start, t_end, ... 
    e, a_max, inertia, w0, q0)
% Junette Hsin 
% Discrete attitude determination for gyrostat 
% 
% Inputs: 
%   dt          = time step 
%   t_start     = time start 
%   t_end       = time end 
%   e           = eigenaxis in inertial frame 
%   a_max       = acceleration constraint 
%   v_max       = velocity constraint 
%   J           = SC inertia, in body frame 
%   w0          = initial velocity 
%   q0          = initial attitude 
% 
% Outputs: 
%   w_out       = output velocity 
%   q_out       = output attitude 
%   torque_out  = output torque 
%   time_out    = output time 

% Ensure eigenaxis is column not row 
if isrow(e) == 1
    e = e'; 
end 

% eigen skew matrix 
e_skew = [  0      -e(3)     e(2);
            e(3)    0       -e(1);
           -e(2)    e(1)     0 ]; 


% Initializing 
nsteps = 10;                % Resolution for attitude solving 
dphi = 0.5*a_max*dt^2;      % FIRST dphi --> how much angle to slew 
index = 1;                  % outputs index 
w = w0; 
q = q0; 
w_out = [ w0'; zeros(length(dt : dt : (t_end - t_start)), 3 )]; 
q_out = [ q0'; zeros(length(dt : dt : (t_end - t_start)), 4 )];
torque_out = [torqueN'; zeros(length(dt : dt : (t_end - t_start)), 3)]; 
time_out = [t_start; zeros(length(dt : dt : (t_end - t_start)), 1)]; 

% Outer for loop, t_start + dt --> t_end 
for t = t_start+dt : dt : t_end 
       
    % rotation around eigenaxis into DCM 
    B_DCM_N = [ cos(dphi)*eye(3) + (1 - cos(dphi))*(e*e') - sin(dphi)*e_skew ]; 
    N_DCM_B = B_DCM_N'; 
    
    % Desired dw from inertial to body 
    dw = B_DCM_N*a_max*e; 
    
    % Solve for torque 
    torque = inertia*dw + cross(w, inertia*w); 
    
    % Solve dynamics and kinematics equations 
    for i = 1:nsteps
        dw = inv(inertia)*(torque - cross(w, inertia*w));
        w_skew = [  0      -w(3)    w(2); 
                    w(3)    0      -w(1); 
                   -w(2)    w(1)    0 ] ; 
        dw = inv(inertia) * ( -w_skew * inertia * w + torque); 
        q_skew = [ q(4)     -q(3)       q(2);
                   q(3)      q(4)      -q(1);
                  -q(2)      q(1)       q(4);
                  -q(1)     -q(2)      -q(3)]; 
        dq = 1/2 * q_skew * w ;
        w = w + dw*dt/nsteps;
        q = q + dq*dt/nsteps;
    end
    
    % Next dphi angle to slew 
    dphi = w*dt + 0.5*a_max*dt^2; 
    
    % Outputs 
    index = index + 1; 
    q_out(index, :) = q'; 
    w_out(index, :) = w'; 
    torque_out(index, :) = torque; 
    time_out(index, 1) = t; 
    
end 
    



