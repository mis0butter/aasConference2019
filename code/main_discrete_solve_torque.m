% MAIN Example 
% Junette Hsin 
% Mohammad Ayoubi 

% clear; 
% close all; 
% main_inputs_P3_G0         % Creates all inputs and variables in workspace 

% optional plotting routine to check things 
plot_option = 0; 
if plot_option == 1
    phi_check_vectors 
end 

%% Determine Phi 1 slew times

% Initial conditions 
w0 = 0; 
t0 = 0; 
wf = 0; 

[t1, t2, t3] = find_slew_times(t0, w0, wf, wMax, aMax, phi1, phi_tt); 

% Discretize slew times to 2 decimal places 
t1 = round(t1, 2); 
t2 = round(t2, 2); 
t3 = round(t3, 3); 

%% Solve for attitude determination - first slew 
        
dt = 1/100; 

% t0 --> t1 
w0 = [    0;    0;      0];                         % wrt G0 frame 
q0 = [    0;    0;      0;      1];                 % wrt G0 frame 
a_P = [ 0;  0; aMax];                               % acceleration around eigenaxis of P frame 
a_G0 = G0_DCM_P*a_P; 
torque_G0 = inertia_SC*a_G0;                                 % wrt G0 frame 

% [t1_phi1, q1_phi1, w1_phi1, torque1_phi1] = gyrostat_discrete(dt, t0, t1, inertia_SC, torque_G0, w0, q0); 
[t1_phi1, q1_phi1, w1_phi1, torque1_phi1] = gyrostat_discrete_torqueN(dt, t0, t1, inertia_SC, torque_G0, w0, q0);

% t1 --> t2 
w0 = w1_phi1(end, :)'; 
q0 = q1_phi1(end, :)'; 
torque = [0; 0; 0]; 

% [t2_phi1, q2_phi1, w2_phi1, torque2_phi1] = gyrostat_discrete(dt, t1, t2, inertia_SC, torque, w0, q0); 
[t2_phi1, q2_phi1, w2_phi1, torque2_phi1] = gyrostat_discrete_torqueN(dt, t1, t2, inertia_SC, torque, w0, q0);

% t2 --> t3 
w0 = w2_phi1(end, :)'; 
q0 = q2_phi1(end, :)'; 

% [t3_phi1, q3_phi1, w3_phi1, torque3_phi1] = gyrostat_discrete(dt, t2, t3, inertia_SC, -torque_G0, w0, q0); 
[t3_phi1, q3_phi1, w3_phi1, torque3_phi1] = gyrostat_discrete_torqueN(dt, t2, t3, inertia_SC, -torque_G0, w0, q0);


t_phi1 = [t1_phi1; t2_phi1(2:end); t3_phi1(2:end)]; 
w_phi1 = [w1_phi1; w2_phi1(2:end ,:); w3_phi1(2:end, :)]; 
q_phi1 = [q1_phi1; q2_phi1(2:end ,:); q3_phi1(2:end, :)]; 
torque_phi1 = [torque1_phi1; torque2_phi1(2:end ,:); torque3_phi1(2:end, :)]; 

%% Plot phi1
plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_phi1, q_phi1, w_phi1, ypr_phi1, 1)

    time_torque = linspace(t0, t3, length(torque_phi1)); 
    figure
        plot(time_torque, torque_phi1); 
        grid on 
        ylabel('Torque (Nm)') 
        xlabel('Time (sec)') 
        title('Torque (Nm)') 
end 
    
%% Determine Phi 2 slew times 

% Starting from rest, ending at rest 
w0 = 0; 
t0 = 0; 
wf = 0; 

% phi2 = 2.5; 
% phi2 = 1.9043;      % hard-coded just for this 
% phi2 = 2.25; 

[t1, t2, t3] = find_slew_times(t0, w0, wf, wMax, aMax, phi2, phi_tt); 

%% Direction of phi2

sun_e = cross(S_G0, e_G0); 
a = sun_e/norm(sun_e); 
b = Pf_G0/norm(Pf_G0); 
c = Pi_G0/norm(Pi_G0); 

% if angle from sun_e and Pf is larger than sun_e and Pi 
if acos(dot(a,b)) > acos(dot(a, c))
    sign = -1; 
else 
    sign = 1; 
end 
    
%% Solve for attitude determination - second slew 
% 
% What I need to happen: 
% for t0 --> t1: 
% 
% the direction of torque_G needs to be recalculated at every time step. 

dt = 1/100; 

% This is the quaternion of current G in initial G0 frame. G0_q_G
q0 = q_phi1(end, :)';               % G0_q_G
w0 = w_phi1(end, :)'; 
phi_w0 = 0; 

[t1_phi2, q1_phi2, w1_phi2, torque1_phi2, phi1_phi2] = gyrostat_discrete_torqueN_solve_torque(dt, t0, t1, ... 
    S_G0, sign*aMax, inertia_SC, w0, q0, phi_w0); 

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

t_phi2 = [t1_phi2; t2_phi2(2:end); t3_phi2(2:end)]; 
w_phi2 = [w1_phi2; w2_phi2(2:end ,:); w3_phi2(2:end, :)]; 
q_phi2 = [q1_phi2; q2_phi2(2:end ,:); q3_phi2(2:end, :)]; 
torque_phi2 = [torque1_phi2; torque2_phi2(2:end, :); torque3_phi2(2:end, :)]; 

%% Plot phi2
plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_phi2, q_phi2, w_phi2, ypr_phi2, 2)
end 

%% Determine Phi 3 slew times

w0 = 0; 
t0 = 0; 
wf = 0; 

[t1, t2, t3] = find_slew_times(t0, w0, wf, wMax, aMax, phi3, phi_tt); 

%% Solve for attitue determination - third slew 
        
dt = 1/100; 

% t0 --> t1 
% tEnd = t1 - t0; 
w0 = w_phi2(end, :)'; 
q0 = q_phi2(end, :)';  
a_P = [0; 0; -aMax]; 
a_G0 = G0_DCM_P*a_P; 
G0_DCM_G = quat2DCM(q0); 
G_DCM_G0 = G0_DCM_G'; 
a_G = G_DCM_G0*a_G0; 
torque_G = inertia_SC*a_G;              % inertia_SC always in G frame 
torque_G0 = G0_DCM_G*torque_G; 

% [t1_phi3, y1_phi3] = ode45(@(t,Z) gyrostat_cont(inertia_SC, torque, Z), [0, tEnd], [w0; q0]);
% [t1_phi3, q1_phi3, w1_phi3, torque1_phi3] = gyrostat_discrete(dt, t0, t1, inertia_SC, torque_G0, w0, q0); 
% [t1_phi3, q1_phi3, w1_phi3, torque1_phi3] = gyrostat_discrete_torqueN(dt, t0, t1, inertia_SC, torque_G0, w0, q0); 
[t1_phi3, q1_phi3, w1_phi3, torque1_phi3, phi1_phi3] = gyrostat_discrete_torqueN_solve_torque(dt, t0, t1, ... 
    e_G0, aMax, inertia_SC, w0, q0, wMax); 

% t1 --> t2 
% tEnd = t2 - t1; 
w0 = w1_phi3(end, :)'; 
q0 = q1_phi3(end, :)'; 
torque = [0; 0; 0]; 

% [t2_phi3, y2_phi3] = ode45(@(t,Z) gyrostat_cont(inertia_SC, torque, Z), [0, tEnd], [w0; q0]); 
% [t2_phi3, q2_phi3, w2_phi3, torque2_phi3] = gyrostat_discrete(dt, t1, t2, inertia_SC, torque, w0, q0); 
% [t2_phi3, q2_phi3, w2_phi3, torque2_phi3] = gyrostat_discrete_torqueN(dt, t1, t2, inertia_SC, torque, w0, q0); 
[t2_phi3, q2_phi3, w2_phi3, torque2_phi3, phi2_phi3] = gyrostat_discrete_torqueN_solve_torque(dt, t1, t2, ... 
    e_G0, 0, inertia_SC, w0, q0, wMax); 

% t2 --> t3 
% tEnd = t3 - t2; 
w0 = w2_phi3(end, :)'; 
q0 = q2_phi3(end, :)'; 

% [t3_phi3, y3_phi3] = ode45(@(t, Z) gyrostat_cont(inertia_SC, torque, Z), [0, tEnd], [w0; q0]); 
% [t3_phi3, q3_phi3, w3_phi3, torque3_phi3] = gyrostat_discrete(dt, t2, t3, inertia_SC, -torque_G, w0, q0); 
% [t3_phi3, q3_phi3, w3_phi3, torque3_phi3] = gyrostat_discrete_torqueN(dt, t2, t3, inertia_SC, -torque_G0, w0, q0); 
[t3_phi3, q3_phi3, w3_phi3, torque3_phi3, phi3_phi3] = gyrostat_discrete_torqueN_solve_torque(dt, t2, t3, ... 
    e_G0, -aMax, inertia_SC, w0, q0, wMax); 

t_phi3 = [t1_phi3; t2_phi3(2:end); t3_phi3(2:end)]; 
w_phi3 = [w1_phi3; w2_phi3(2:end ,:); w3_phi3(2:end, :)]; 
q_phi3 = [q1_phi3; q2_phi3(2:end ,:); q3_phi3(2:end, :)]; 
torque_phi3 = [torque1_phi3; torque2_phi3(2:end ,:); torque3_phi3(2:end, :)]; 

%% Plot phi3
plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_phi3, q_phi3, w_phi3, ypr_phi3, 3)
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

plot_option = 1; 
if plot_option == 1
    plot_qwypr(t_total, q_total, w_total, torque_total, a_total, 'total', aMax, wMax)
end 

%% PHI NOMINAL - IF THERE WAS NO SUN INTRUSION 

% Determine Phi Nom slew times

% Initial conditions 
w0 = 0; 
t0 = 0; 
wf = 0; 

[t1, t2, t3] = find_slew_times(t0, w0, wf, wMax, aMax, theta_Pi_Pf, phi_tt); 

% Discretize slew times to 2 decimal places 
t1 = round(t1, 2); 
t2 = round(t2, 2); 
t3 = round(t3, 3); 

%% Solve for attitude determination - phiNom
        
% t0 --> t1 
w0 = [    0;    0;      0];                         % wrt G0 frame 
q0 = [    0;    0;      0;      1];                 % wrt G0 frame 
a_P = [ 0;  0; aMax];                               % acceleration around eigenaxis of P frame 
a_G0 = G0_DCM_P*a_P; 
torque_G0 = inertia_SC*a_G0;                                 % wrt G0 frame 

dt = 1/100; 
[t1_phiNom, q1_phiNom, w1_phiNom, torque1_phiNom] = gyrostat_discrete_torqueN(dt, t0, t1, inertia_SC, torque_G0, w0, q0);

% t1 --> t2 
w0 = w1_phiNom(end, :)'; 
q0 = q1_phiNom(end, :)'; 
torque = [0; 0; 0]; 

[t2_phiNom, q2_phiNom, w2_phiNom, torque2_phiNom] = gyrostat_discrete_torqueN(dt, t1, t2, inertia_SC, torque, w0, q0);

% t2 --> t3 
w0 = w2_phiNom(end, :)'; 
q0 = q2_phiNom(end, :)'; 
% torque = -torque_G0; 

[t3_phiNom, q3_phiNom, w3_phiNom, torque3_phiNom] = gyrostat_discrete_torqueN(dt, t2, t3, inertia_SC, -torque_G0, w0, q0);

t_phiNom = [t1_phiNom; t2_phiNom(2:end); t3_phiNom(2:end)]; 
w_phiNom = [w1_phiNom; w2_phiNom(2:end ,:); w3_phiNom(2:end, :)]; 
q_phiNom = [q1_phiNom; q2_phiNom(2:end ,:); q3_phiNom(2:end, :)]; 
torque_phiNom = [torque1_phiNom; torque2_phiNom(2:end ,:); torque3_phiNom(2:end, :)]; 

%% Post processing 

post_processing_N

