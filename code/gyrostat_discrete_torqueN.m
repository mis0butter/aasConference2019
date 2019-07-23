function [time, q_out, w_out, torque_out] = gyrostat_discrete_torqueN(dt, t_start, t_end, inertia, torqueN, w0, q0)
% Junette Hsin 
% Discrete attitude determination for gyrostat 
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

% w and Q need to be columns 
if isrow(w0) == 1
    w0 = w0'; 
end 
if isrow(q0) == 1
    q0 = q0'; 
end 
    
% Initialize 
w = w0; 
q = q0;                     % G0_q_G

N_DCM_G = quat2DCM(q);          % G0_DCM_G
torque = N_DCM_G'*torqueN; 

int = 1; 
w_out = [ w0'; zeros(length(dt : dt : (t_end - t_start)), 3 )]; 
q_out = [ q0'; zeros(length(dt : dt : (t_end - t_start)), 4 )];
torque_out = [torqueN'; zeros(length(dt : dt : (t_end - t_start)), 3)]; 
time = [t_start; zeros(length(dt : dt : (t_end - t_start)), 1)]; 

% Outer for loop, t_start + dt --> t_end 
for t = t_start+dt : dt : t_end 
% for t = t_start : dt : t_end 

    if norm(torque) == 0

        % do nothing; don't calculate dw or dq or new torque 
        
    else 

        nsteps = 10;
        q_old = q;      % G0_q_G(old)

        % Break into steps for attitude determination 
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
    %         w = w + dw;
    %         q = q + dq;
        end
        
        N_DCM_G = quat2DCM(q);
        torque = N_DCM_G'*torqueN; 
        
    end 
    
    int = int + 1; 
    q_out(int, :) = q'; 
    w_out(int, :) = w'; 
    torque_out(int, :) = torque; 
    time(int) = t; 
    
end 
end 

% inputs: dt, inertia, torque, w0, q0
% outputs: time, q, w 
