function [time, q_out, w_out] = gyrostat_discrete(dt, t_start, t_end, inertia, torque, w0, q0)

% w and Q need to be columns 
if isrow(w0) == 1
    w0 = w0'; 
end 
if isrow(q0) == 1
    q0 = q0'; 
end 
    
w = w0; 
q = q0; 
int = 1; 
w_out = [ w0'; zeros(length(dt : dt : (t_end - t_start)), 3 )]; 
q_out = [ q0'; zeros(length(dt : dt : (t_end - t_start)), 4 )];
time = [t_start; zeros(length(dt : dt : (t_end - t_start)), 1)]; 

for t = t_start+dt : dt : t_end 
    
    nsteps = 10;
    
    for i = 1:nsteps
        dw = inv(inertia)*(torque - cross(w, inertia*w));
        dq = 1/2*[q(4), -q(3), q(2);
                  q(3), q(4), -q(1);
                  -q(2), q(1), q(4);
                  -q(1), -q(2), -q(3)]*w;
        w = w + dw*dt/nsteps;
        q = q + dq*dt/nsteps;
    end
    
    DCM = SpinCalc('QtoDCM', q', eps, 0); 
    torque = DCM*torque; 
%     inertia = DCM*inertia; 
    
    int = int + 1; 
    q_out(int, :) = q'; 
    w_out(int, :) = w'; 
    time(int) = t; 
    
end 

% inputs: dt, inertia, torque, w0, q0
% outputs: time, q, w 

end 