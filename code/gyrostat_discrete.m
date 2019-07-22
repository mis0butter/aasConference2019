function [time, q_out, w_out] = gyrostat_discrete(dt, t_start, t_end, inertia, torque, w0, q0)
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
q = q0; 
int = 1; 
w_out = [ w0'; zeros(length(dt : dt : (t_end - t_start)), 3 )]; 
q_out = [ q0'; zeros(length(dt : dt : (t_end - t_start)), 4 )];
time = [t_start; zeros(length(dt : dt : (t_end - t_start)), 1)]; 

% Outer for loop, t_start + dt --> t_end 
for t = t_start+dt : dt : t_end 
% for t = t_start : dt : t_end 
    
    nsteps = 10;
    
    % Break into steps for attitude determination 
    for i = 1:nsteps
%         dw = inv(inertia)*(torque - cross(w, inertia*w));
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
    
    % Find rotation from previous 
    DCM = SpinCalc('QtoDCM', q', eps, 0); 
%     torque = DCM*torque; 
%     inertia = DCM*inertia; 
    
    int = int + 1; 
    q_out(int, :) = q'; 
    w_out(int, :) = w'; 
    time(int) = t; 
    
end 

% inputs: dt, inertia, torque, w0, q0
% outputs: time, q, w 

end 