function [t_out, q_out, w_out, torque_out] = gyrostat_CD(dt, t_start, t_end, inertia, torque0, w0, q0)
% Junette Hsin 
% Continuous-Discrete attitude determination for gyrostat 
% 
% Inputs: 
%   dt          = time step (s) 
%   t_start     = start time in slew 
%   t_end       = end time in slew 
%   inertia     = inertia in gyrostat frame 
%   torque      = torque in gyrostat frame 
%   w0          = initial angular velocity 
%   q0          = initial quaternions 
% 
% Outputs: 
%   time        = time vector at dt steps 
%   q_out       = output quaternions
%   w_out       = output angular velocity 
%   torque_out  = capture torque 

% Initial inputs 
w_in = w0; 
q_in = q0; 
torque = torque0; 

% Initialize outputs
t_out = 0; 
w_out = w_in'; 
q_out = q_in'; 
torque_out = torque; 

% for t_i = 0 : dt : tEnd - dt 
for t_i = t_start+dt : dt : t_end
    
    q_old = q_in; 
    
    [t, y] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [0, dt], [w_in; q_in]);
    w = y(:, 1:3); 
    q = y(:, 4:7); 
    w_in = w(end, :)'; 
    q_in = q(end, :)'; 
    
    q_new = q_in;      % G0_q_G(new) 
    old_q_new = quat_diff(q_old, q_new);    % B_q_C = quat_diff(A_q_B, A_q_C)   
    old_DCM_new = quat2DCM(old_q_new); 
    new_DCM_old = old_DCM_new'; 
    torque = new_DCM_old*torque; 
    
    if isrow(t) == 1
        t = t'; 
    end 
    t_out = [t_out; t_out(end) + t(2:end)]; 
    w_out = [w_out; w(2:end, :)]; 
    q_out = [q_out; q(2:end, :)]; 
    torque_out = [torque_out; torque(2:end, :)]; 

end