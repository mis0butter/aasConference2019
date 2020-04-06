function [time_out, q_out, w_out, h_out, torque_out, ws_out, phi] = propagate_attitude ... 
    (dt, t_start, t_end, e, a_max, I_sc, w, q, h, w_max)
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

%% hacks 
% wheel inertia 
I_whl = 0.0019*(eye(4)); 

% alpha = 60 deg 
a = 65*pi/180; 

% Wheel axes from SC axes transformation matrix [4x3] 
whl_mx_sc = [ ... 
    sqrt(2)/2 * sin(a)      -sqrt(2)/2 * sin(a)     cos(a)  ; 
    sqrt(2)/2 * sin(a)      sqrt(2)/2 * sin(a)      cos(a)  ; 
    -sqrt(2)/2 * sin(a)     -sqrt(2)/2 * sin(a)     cos(a)  ; 
    -sqrt(2)/2 * sin(a)     sqrt(2)/2 * sin(a)      cos(a)  ]; 

% SC axes from wheel transformation matrix [3x4] 
sc_mx_whl = (whl_mx_sc' * whl_mx_sc) ^-1 * whl_mx_sc'; 
A = sc_mx_whl;  % for brevity 

%%

% Initializing 
nsteps = 10;                % Resolution for attitude solving 
phi = 0;                    % FIRST dphi --> how much angle to slew 
index = 1;                  % outputs index 
w_out = [ w'; zeros(length(dt : dt : (t_end - t_start)), 3 )]; 
q_out = [ q'; zeros(length(dt : dt : (t_end - t_start)), 4 )];
h_out = [ h'; zeros(length(dt : dt : (t_end - t_start)), 3 )];
torque_out = zeros(length(dt : dt : (t_end - t_start)) + 1, 3); 
time_out = [t_start; zeros(length(dt : dt : (t_end - t_start)), 1)]; 
ws_out = zeros(size(q_out)); 

% Outer for loop, t_start + dt --> t_end 
for t = t_start+dt : dt : t_end 
    
    % Notes from meeting with Daniel 
    % torque limit from wheel - assume max torque per wheel 
    % calculate a_max from current SC orientation. how much torque can you
    % get from each direction? 
    
    % Desired phi angle 
    phi = w_max*(t - t_start) + 0.5*a_max*(t - t_start)^2; 
       
    % rotation around eigenaxis by phi into DCM 
    B_DCM_N = [ cos(phi).*eye(3) + (1 - cos(phi)).*(e*e') - sin(phi).*e_skew ]; 
    N_DCM_B = B_DCM_N'; 
    
    % Desired dw from inertial to body 
    dw = B_DCM_N*a_max*e; 
    
    % Solve for torque - make sure it does not exceed max torque in that
    % direction 
    torque = I_sc*dw + cross(w, I_sc*w); 
    
%     if torque exceeds maximum allowed torque by the reaction wheels 
%         decrease this torque vector so that it is < threshold, keep the direction 
%     end 
    
    % Solve dynamics and kinematics equations 
    for i = 1:nsteps
%         dw = inv(inertia)*(torque - cross(w, inertia*w));
        w_skew = [  0      -w(3)    w(2); 
                    w(3)    0      -w(1); 
                   -w(2)    w(1)    0 ] ; 
        q_skew = [ q(4)     -q(3)       q(2);
                   q(3)      q(4)      -q(1);
                  -q(2)      q(1)       q(4);
                  -q(1)     -q(2)      -q(3)]; 
        dw = inv(I_sc) * ( -w_skew * I_sc * w + torque); 
        dq = 1/2 * q_skew * w ;
        dh = -torque - w_skew * h; 
        w = w + dw*dt/nsteps;
        q = q + dq*dt/nsteps;
        h = h + dh*dt/nsteps; 
        
        % Wheel speeds 
        ws = inv(I_whl) * ( ... 
            (A'*A)^-1 * A' * h - I_whl * A' * w ... 
            ); 
    end 
    
    % Outputs 
    index = index + 1; 
    q_out(index, :) = q'; 
    w_out(index, :) = w'; 
    h_out(index, :) = h'; 
    torque_out(index, :) = torque; 
    time_out(index, 1) = t; 
    ws_out(index, :) = ws; 
    
end 
    



