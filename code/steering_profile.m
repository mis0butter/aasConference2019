% STEERING PROFILE
% Junette Hsin

dt = 1/100;                 % gyrostat attitude solver discretization 

%% Determine Phi 1 slew times

% Initial conditions 
w0 = 0; 
t0 = 0; 
wf = 0; 

[t1, t2, t3] = find_slew_times(t0, w0, wf, wMax, aMax, phi1, phi_t); 

% Discretize slew times to 2 decimal places 
t1 = round(t1, 2); 
t2 = round(t2, 2); 
t3 = round(t3, 3); 

% %% Solve for attitude determination - first slew 
% 
% % t0 --> t1 
% w0 = [    0;    0;      0];                         % wrt G0 frame 
% q0 = [    0;    0;      0;      1];                 % wrt G0 frame 
% a_P = [ 0;  0; aMax];                               % acceleration around eigenaxis of P frame 
% a_G0 = G0_DCM_P*a_P; 
% torque_G0 = inertia_SC*a_G0;                                 % wrt G0 frame 
% 
% [t1_phi1, q1_phi1, w1_phi1, torque1_phi1] = gyrostat_discrete_torqueN(dt, t0, t1, inertia_SC, torque_G0, w0, q0);
% 
% % t1 --> t2 
% w0 = w1_phi1(end, :)'; 
% q0 = q1_phi1(end, :)'; 
% torque = [0; 0; 0]; 
% 
% [t2_phi1, q2_phi1, w2_phi1, torque2_phi1] = gyrostat_discrete_torqueN(dt, t1, t2, inertia_SC, torque, w0, q0);
% 
% % t2 --> t3 
% w0 = w2_phi1(end, :)'; 
% q0 = q2_phi1(end, :)'; 
% 
% [t3_phi1, q3_phi1, w3_phi1, torque3_phi1] = gyrostat_discrete_torqueN(dt, t2, t3, inertia_SC, -torque_G0, w0, q0);
% 
% % Putting it all together 
% t_phi1 = [t1_phi1; t2_phi1(2:end); t3_phi1(2:end)]; 
% w_phi1 = [w1_phi1; w2_phi1(2:end ,:); w3_phi1(2:end, :)]; 
% q_phi1 = [q1_phi1; q2_phi1(2:end ,:); q3_phi1(2:end, :)]; 
% torque_phi1 = [torque1_phi1; torque2_phi1(2:end ,:); torque3_phi1(2:end, :)]; 

%% Solve for attitude determination - first slew 

% t0 --> t1 
w0 = [    0;    0;      0];                         % wrt G0 frame 
q0 = [    0;    0;      0;      1];                 % wrt G0 frame 

[t1_phi1, q1_phi1, w1_phi1, torque1_phi1, phi1_phi1] = gyrostat_discrete_torqueN_solve_torque(dt, t0, t1, ... 
    e_G0, aMax, inertia_SC, w0, q0, wMax); 

% t1 --> t2 
w0 = w1_phi1(end, :)'; 
q0 = q1_phi1(end, :)'; 

[t2_phi1, q2_phi1, w2_phi1, torque2_phi1, phi2_phi1] = gyrostat_discrete_torqueN_solve_torque(dt, t1, t2, ... 
    e_G0, 0, inertia_SC, w0, q0, wMax); 

% t2 --> t3 
w0 = w2_phi1(end, :)'; 
q0 = q2_phi1(end, :)'; 

[t3_phi1, q3_phi1, w3_phi1, torque3_phi1, phi3_phi1] = gyrostat_discrete_torqueN_solve_torque(dt, t2, t3, ... 
    e_G0, -aMax, inertia_SC, w0, q0, wMax); 

% Putting it all together 
t_phi1 = [t1_phi1; t2_phi1(2:end); t3_phi1(2:end)]; 
w_phi1 = [w1_phi1; w2_phi1(2:end ,:); w3_phi1(2:end, :)]; 
q_phi1 = [q1_phi1; q2_phi1(2:end ,:); q3_phi1(2:end, :)]; 
torque_phi1 = [torque1_phi1; torque2_phi1(2:end ,:); torque3_phi1(2:end, :)]; 

%% Plot phi1

% acceleration stuff 
a_phi1 = zeros(length(w_phi1) - 1, 3); 
for i = 1:length(w_phi1) - 1 
    a_phi1(i, :) = (1/dt)*(w_phi1(i + 1, :) - w_phi1(i, :)); 
end 

plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_phi1, q_phi1, w_phi1, torque_phi1, a_phi1, 'total', aMax, wMax)
end 
    
%% Determine Phi 2 slew times 

% Starting from rest, ending at rest 
w0 = 0; 
t0 = 0; 
wf = 0; 

[t1, t2, t3] = find_slew_times(t0, w0, wf, wMax, aMax, phi2, phi_t); 

%% Direction of phi2

% depends on angle between sun and eigenaxis 
if acos(dot(e_G0, S_G0)) < pi/2
    sign = 1; 
else 
    sign = -1; 
end 
    
%% Solve for attitude determination - second slew 

% This is the quaternion of current G in initial G0 frame. G0_q_G
q0 = q_phi1(end, :)';               % G0_q_G
w0 = w_phi1(end, :)'; 
phi_w0 = 0; 

[t1_phi2, q1_phi2, w1_phi2, torque1_phi2, phi1_phi2] = gyrostat_discrete_torqueN_solve_torque(dt, t0, t1, ... 
    S_G0, sign*aMax, inertia_SC, w0, q0, wMax); 

% t1 --> t2 
w0 = w1_phi2(end, :)'; 
q0 = q1_phi2(end, :)'; 

[t2_phi2, q2_phi2, w2_phi2, torque2_phi2, phi2_phi2] = gyrostat_discrete_torqueN_solve_torque(dt, t1, t2, ... 
    S_G0, 0, inertia_SC, w0, q0, wMax); 

% t2 --> t3 
w0 = w2_phi2(end, :)'; 
q0 = q2_phi2(end, :)';  

[t3_phi2, q3_phi2, w3_phi2, torque3_phi2, phi3_phi2] = gyrostat_discrete_torqueN_solve_torque(dt, t2, t3, ... 
    S_G0, -sign*aMax, inertia_SC, w0, q0, wMax); 

% Putting it all together 
t_phi2 = [t1_phi2; t2_phi2(2:end); t3_phi2(2:end)]; 
w_phi2 = [w1_phi2; w2_phi2(2:end ,:); w3_phi2(2:end, :)]; 
q_phi2 = [q1_phi2; q2_phi2(2:end ,:); q3_phi2(2:end, :)]; 
torque_phi2 = [torque1_phi2; torque2_phi2(2:end, :); torque3_phi2(2:end, :)]; 

%% Plot phi2

% acceleration stuff 
a_phi2 = zeros(length(w_phi2) - 1, 3); 
for i = 1:length(w_phi2) - 1 
    a_phi2(i, :) = (1/dt)*(w_phi2(i + 1, :) - w_phi2(i, :)); 
end 

plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_phi2, q_phi2, w_phi2, torque_phi2, a_phi2, 'total', aMax, wMax)
end 

%% Determine Phi 3 slew times

w0 = 0; 
t0 = 0; 
wf = 0; 

[t1, t2, t3] = find_slew_times(t0, w0, wf, wMax, aMax, phi3, phi_t); 

%% Solve for attitue determination - third slew 

% t0 --> t1 
w0 = w_phi2(end, :)'; 
q0 = q_phi2(end, :)';  

[t1_phi3, q1_phi3, w1_phi3, torque1_phi3, phi1_phi3] = gyrostat_discrete_torqueN_solve_torque(dt, t0, t1, ... 
    e_G0, aMax, inertia_SC, w0, q0, wMax); 

% t1 --> t2 
w0 = w1_phi3(end, :)'; 
q0 = q1_phi3(end, :)'; 

[t2_phi3, q2_phi3, w2_phi3, torque2_phi3, phi2_phi3] = gyrostat_discrete_torqueN_solve_torque(dt, t1, t2, ... 
    e_G0, 0, inertia_SC, w0, q0, wMax); 

% t2 --> t3 
w0 = w2_phi3(end, :)'; 
q0 = q2_phi3(end, :)'; 

[t3_phi3, q3_phi3, w3_phi3, torque3_phi3, phi3_phi3] = gyrostat_discrete_torqueN_solve_torque(dt, t2, t3, ... 
    e_G0, -aMax, inertia_SC, w0, q0, wMax); 

% Putting it all together 
t_phi3 = [t1_phi3; t2_phi3(2:end); t3_phi3(2:end)]; 
w_phi3 = [w1_phi3; w2_phi3(2:end ,:); w3_phi3(2:end, :)]; 
q_phi3 = [q1_phi3; q2_phi3(2:end ,:); q3_phi3(2:end, :)]; 
torque_phi3 = [torque1_phi3; torque2_phi3(2:end ,:); torque3_phi3(2:end, :)]; 

%% Plot phi3

% acceleration stuff 
a_phi3 = zeros(length(w_phi3) - 1, 3); 
for i = 1:length(w_phi3) - 1 
    a_phi3(i, :) = (1/dt)*(w_phi3(i + 1, :) - w_phi3(i, :)); 
end 

plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_phi3, q_phi3, w_phi3, torque_phi3, a_phi3, 'total', aMax, wMax)
end 

%% total slew stuff 

t_total = [ t_phi1; ... 
            t_phi1(end) + t_phi2(2:end); ... 
            t_phi1(end) + t_phi2(end) + t_phi3(2:end)]; 
        
w_total = [ w_phi1; ... 
            w_phi2(2:end, :); ... 
            w_phi3(2:end, :)]; 
        
q_total = [ q_phi1; ... 
            q_phi2(2:end, :); ...  
            q_phi3(2:end, :)]; 
        
torque_total = [torque_phi1; ... 
            torque_phi2(2:end, :); ... 
            torque_phi3(2:end, :)]; 
        
%% plot total stuff 

% acceleration stuff 
a_total = zeros(length(w_total) - 1, 3); 
for i = 1:length(w_total) - 1 
    a_total(i, :) = (1/dt)*(w_total(i + 1, :) - w_total(i, :)); 
end 

plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_total, q_total, w_total, torque_total, a_total, 'total', aMax, wMax)
end 

%% PHI NOMINAL - IF THERE WAS NO SUN INTRUSION 

% Determine Phi Nom slew times

% Initial conditions 
w0 = 0; 
t0 = 0; 
wf = 0; 

[t1, t2, t3] = find_slew_times(t0, w0, wf, wMax, aMax, theta_Pi_Pf, phi_t); 

% Discretize slew times to 2 decimal places 
t1 = round(t1, 2); 
t2 = round(t2, 2); 
t3 = round(t3, 3); 

%% Solve for attitude determination - phiNom
        
% t0 --> t1 
w0 = [    0;    0;      0];                         % wrt G0 frame 
q0 = [    0;    0;      0;      1];                 % wrt G0 frame 

[t1_phiNom, q1_phiNom, w1_phiNom, torque1_phiNom] = gyrostat_discrete_torqueN(dt, t0, t1, inertia_SC, torque_G0, w0, q0);

% t1 --> t2 
w0 = w1_phiNom(end, :)'; 
q0 = q1_phiNom(end, :)'; 
torque = [0; 0; 0]; 

[t2_phiNom, q2_phiNom, w2_phiNom, torque2_phiNom] = gyrostat_discrete_torqueN(dt, t1, t2, inertia_SC, torque, w0, q0);

% t2 --> t3 
w0 = w2_phiNom(end, :)'; 
q0 = q2_phiNom(end, :)'; 

[t3_phiNom, q3_phiNom, w3_phiNom, torque3_phiNom] = gyrostat_discrete_torqueN(dt, t2, t3, inertia_SC, -torque_G0, w0, q0);

t_phiNom = [t1_phiNom; t2_phiNom(2:end); t3_phiNom(2:end)]; 
w_phiNom = [w1_phiNom; w2_phiNom(2:end ,:); w3_phiNom(2:end, :)]; 
q_phiNom = [q1_phiNom; q2_phiNom(2:end ,:); q3_phiNom(2:end, :)]; 
torque_phiNom = [torque1_phiNom; torque2_phiNom(2:end ,:); torque3_phiNom(2:end, :)]; 