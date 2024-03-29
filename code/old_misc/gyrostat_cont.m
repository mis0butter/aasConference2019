function dZdt = gyrostat_cont(inertia, torque, Z)
% Attitude determination function for finding w and q from inertia, torque,
% and initial w and q 
% 
% Inputs: 
% inertia = [3x3] 
% torque = [3x1] 
% Z(1 - 3) = w
% Z(4 - 7) = q

% Ensure torque is column [3x1]
if isrow(torque) == 1
    torque = torque'; 
end 
    
w = [Z(1); Z(2); Z(3)]; 
q = [Z(4); Z(5); Z(6); Z(7)]; 

w_skew = [  0      -w(3)    w(2); 
            w(3)    0      -w(1); 
           -w(2)    w(1)    0 ] ; 
dw = inv(inertia) * ( -w_skew * inertia * w + torque); 
q_skew = [ q(4)     -q(3)       q(2);
           q(3)      q(4)      -q(1);
          -q(2)      q(1)       q(4);
          -q(1)     -q(2)      -q(3)]; 
dq = 1/2 * q_skew * w ;

dZdt = [dw; dq]; 

% Previous formula 
% 
% w_skew =  [ 0      -Z(3)   Z(2); 
%             Z(3)    0     -Z(1); 
%            -Z(2)    Z(1)   0]; 
%         
% q_skew = [  0       Z(3)   -Z(2)    Z(1); 
%            -Z(3)    0      -Z(1)    Z(2); 
%             Z(2)   -Z(1)    0       Z(3); 
%            -Z(1)    Z(2)   -Z(3)    0]; 
%        
% % w dot [3x1] = ... 
% % q dot [4x1] = ... 
% 
% dZdt = [    inv(inertia)*(-w_skew * inertia * [Z(1); Z(2); Z(3)] + torque) ; 
%             0.5 * q_skew * [Z(4); Z(5); Z(6); Z(7)] ]; 
    


% dw = inv(inertia)*( ... 
%     - [ 0      -w(3)    w(2); 
%         w(3)    0      -w(1); 
%        -w(2)    w(1)    0 ] ... 
%     * inertia * [Z(1); Z(2); Z(3)] + torque);
% dq = 1/2* ... 
%     [ q(4)     -q(3)       q(2);
%       q(3)      q(4)      -q(1);
%      -q(2)      q(1)       q(4);
%      -q(1)     -q(2)      -q(3)] * w ;
% 
% dq = 0.5*[  0       w(3)   -w(2)    w(1);     % Interesting - this 
%            -w(3)    0      -w(1)    w(2);     % skew matrix results in 
%             w(2)   -w(1)    0       w(3);     % q's blowing up 
%            -w(1)    w(2)   -w(3)    0]*q; 
        