function dZdt = gyrostat_cont(inertia, torque, Z)

% Z(1 - 3) = w
% Z(4 - 7) = q

% inertia = [408 0 0; 0 427 0; 0 0 305]; 

w_skew =  [ 0      -Z(3)   Z(2); 
            Z(3)    0     -Z(1); 
           -Z(2)    Z(1)   0]; 
        
q_skew = [  0       Z(3)   -Z(2)    Z(1); 
           -Z(3)    0      -Z(1)    Z(2); 
            Z(2)   -Z(1)    0       Z(3); 
           -Z(1)    Z(2)   -Z(3)    0]; 
       
% w dot [3x1] = ... 
% q dot [4x1] = ... 

dZdt = [    inv(inertia)*(-w_skew * inertia * [Z(1); Z(2); Z(3)] + torque) ; 
            0.5 * q_skew * [Z(4); Z(5); Z(6); Z(7)] ]; 

        