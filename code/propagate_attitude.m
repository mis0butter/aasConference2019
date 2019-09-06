function [time_out, q_out, w_out, torque_out, phi] = propagate_attitude(dt, t_start, t_end, ... 
    e, a_max, inertia, w0, q0, w_max)
% 
% Junette Hsin 
% Discrete attitude determination & propagation function for gyrostat 
% 
% Inputs: 
%   dt          = time step 
%   t_start     = time start 
%   t_end       = time end 
%   e           = eigenaxis (axis of rotation) in inertial frame 
%   a_max       = acceleration constraint 
%   inertia     = SC inertia, in body frame 
%   w0          = initial velocity 
%   q0          = initial attitude 
%   w_max       = velocity constraint 
% 
% Outputs: 
%   time_out    = output time 
%   q_out       = output attitude 
%   w_out       = output velocity 
%   torque_out  = output torque 
%   phi         = output phi angle 

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
phi = 0;                    % FIRST dphi --> how much angle to slew 
index = 1;                  % outputs index 
w = w0; 
q = q0; 
w_out = [ w0'; zeros(length(dt : dt : (t_end - t_start)), 3 )]; 
q_out = [ q0'; zeros(length(dt : dt : (t_end - t_start)), 4 )];
torque_out = zeros(length(dt : dt : (t_end - t_start)) + 1, 3); 
time_out = [t_start; zeros(length(dt : dt : (t_end - t_start)), 1)]; 

% Outer for loop, t_start + dt --> t_end 
for t = t_start+dt : dt : t_end 
    
    % Desired phi angle 
    phi = w_max*(t - t_start) + 0.5*a_max*(t - t_start)^2; 
       
    % rotation around eigenaxis by phi into DCM 
    B_DCM_N = [ cos(phi).*eye(3) + (1 - cos(phi)).*(e*e') - sin(phi).*e_skew ]; 
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
    
    % Outputs 
    index = index + 1; 
    q_out(index, :) = q'; 
    w_out(index, :) = w'; 
    torque_out(index, :) = torque; 
    time_out(index, 1) = t; 
    
end 
    



