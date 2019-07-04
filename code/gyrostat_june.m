function [w_out, q_out] = gyrostat_june(torque, inertia, w_in, q_in, dt)

% Not the correct equation below: 
% dw = inv(inertia)*torque; 

% Equation below for angular velocity: 
dw = -[ 0           -w_in(3)   w_in(2); 
        w_in(3)     0           -w_in(1); 
        -w_in(2)    w_in(1)     0]*inertia*w_in; 

% 
dQ = 0.5 * [q_in(4)    -q_in(3)   q_in(2); 
            q_in(3)    q_in(4)    -q_in(1);
            -q_in(2)   q_in(1)    q_in(4); 
            -q_in(1)   -q_in(2)   -q_in(3)];
        
w_out = w_in + dw*dt; 

dQ_mult = [ dQ(4)   dQ(3)   -dQ(2)  dQ(1); 
            -dQ(3)  dQ(4)   dQ(1)   dQ(2); 
            dQ(2)   -dQ(1)  dQ(4)   dQ(3); 
            -dQ(1)  -dQ(2)  -dQ(3)  dQ(4)]; 

q_out = dQ_mult*q_in*dt; 

%%



              
              

                  
