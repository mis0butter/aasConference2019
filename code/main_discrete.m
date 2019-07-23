% MAIN Example 
% Junette Hsin 
% Mohammad Ayoubi 

clear; 
close all; 
main_inputs             % Creates all inputs and variables in workspace 

% optional plotting routine to check things 
plot_option = 0; 
phi_check_vectors 

%% Determine Phi 1 slew times

% Initial conditions 
w0 = 0; 
t0 = 0; 
wf = 0; 
[t1, t2, t3] = find_slew_times(aMax, wMax, phi1, w0, wf, t0); 

% Discretize slew times to 2 decimal places 
t1 = round(t1, 2); 
t2 = round(t2, 2); 
t3 = round(t3, 3); 

%% Solve for attitude determination - first slew 
        
% t0 --> t1 
w0 = [    0;    0;      0];                         % G0_w_G
q0 = [    0;    0;      0;      1];                 % G0_q_G
a = [ 0;  0; aMax];  
torque = inertia*a;                                 % wrt G frame 

% dt = 1/100; 
% [t1_phi1, q1_phi1, w1_phi1, torque1_phi1] = gyrostat_discrete(dt, t0, t1, inertia, torque, w0, q0); 
[t1_phi1, y1_phi1] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [0, t1], [w0; q0]);
w1_phi1 = y1_phi1(:, 1:3); 
q1_phi1 = y1_phi1(:, 4:7); 

% t1 --> t2 
w0 = w1_phi1(end, :)'; 
q0 = q1_phi1(end, :)'; 
a = [0; 0; 0]; 
torque = inertia*a; 

% [t2_phi1, q2_phi1, w2_phi1, torque2_phi1] = gyrostat_discrete(dt, t1, t2, inertia, torque, w0, q0); 
[t2_phi1, y2_phi1] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [t1, t2], [w0; q0]);
w2_phi1 = y2_phi1(:, 1:3); 
q2_phi1 = y2_phi1(:, 4:7); 

% t2 --> t3 
w0 = w2_phi1(end, :)'; 
q0 = q2_phi1(end, :)'; 
a = [ 0; 0; -aMax]; 
torque = inertia*a; 

% [t3_phi1, q3_phi1, w3_phi1, torque3_phi1] = gyrostat_discrete(dt, t2, t3+dt, inertia, torque, w0, q0);
[t3_phi1, y3_phi1] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [t2, t3], [w0; q0]);
w3_phi1 = y3_phi1(:, 1:3); 
q3_phi1 = y3_phi1(:, 4:7);  

% t_phi1 = [t1_phi1; t2_phi1(2:end); t3_phi1(2:end)];
t_phi1 = [t1_phi1; t2_phi1(2:end); t3_phi1(2:end)]; 
w_phi1 = [w1_phi1; w2_phi1(2:end ,:); w3_phi1(2:end, :)]; 
q_phi1 = [q1_phi1; q2_phi1(2:end ,:); q3_phi1(2:end, :)]; 
% torque_phi1 = [torque1_phi1; torque2_phi1(2:end ,:); torque3_phi1(2:end, :)]; 

ypr_phi1 = zeros(length(q_phi1), 3); 
for i = 1:max(size(q_phi1))
    ypr_phi1(i, :) = SpinCalc('QtoEA321', q_phi1(i, :), eps, 0); 
end 

%% Plot phi1
plot_option = 1; 
if plot_option == 1
    plot_qwypr(t_phi1, q_phi1, w_phi1, ypr_phi1, 1)

%     time_torque = linspace(t0, t3, length(torque_phi1)); 
%     figure
%         plot(time_torque, torque_phi1); 
%         grid on 
%         ylabel('Torque (Nm)') 
%         xlabel('Time (sec)') 
%         title('Torque (Nm)') 
end 

%%%%%%
% IF angular separation is less than or equal to payload half-cone angle --> find phi2
% and phi3. Otherwise, slew is just phi1. 
%%%%%%
if alpha > ep
    
    t_phi2 = 0; 
    q_phi2 = zeros(1, 4); 
    w_phi2 = zeros(1, 3); 
    
    t_phi3 = 0; 
    q_phi3 = zeros(1, 4); 
    w_phi3 = zeros(1, 3); 
    
else 
    
    %% Determine Phi 2 slew times 

    % Starting from rest, ending at rest 
    w0 = 0; 
    t0 = 0; 
    wf = 0; 
    [t1, t2, t3] = find_slew_times(aMax, wMax, phi2, w0, wf, t0); 

    % Discretize slew times to 2 decimal places 
    t1 = round(t1, 2); 
    t2 = round(t2, 2); 
    t3 = round(t3, 3); 

    %% Solve for attitude determination - second slew 
    % 
    % What I need to happen: 
    % for t0 --> t1: 
    % 
    % the direction of torque_G needs to be recalculated at every time step. 

    dt = 1/100; 

    % This is the quaternion of current G in initial G0 frame. G0_q_G
    q0 = q_phi1(end, :)';    
    torque_G0 = inertia*aMax*S_G0; 

    w0 = w_phi1(end, :)'; 

    [t1_phi2, q1_phi2, w1_phi2, torque1_phi2] = gyrostat_discrete_torqueN(dt, t0, t1, inertia, torque_G0, w0, q0); 

    % t1 --> t2 
    w0 = w1_phi2(end, :)'; 
    q0 = q1_phi2(end, :)'; 
    a = [0; 0; 0]; 
    torque = inertia*a; 

    [t2_phi2, q2_phi2, w2_phi2, torque2_phi2] = gyrostat_discrete_torqueN(dt, t1, t2, inertia, torque, w0, q0); 

    % t2 --> t3 
    w0 = w2_phi2(end, :)'; 
    q0 = q2_phi2(end, :)'; 
    torque_G0 = inertia*-aMax*S_G0; 

    % a = -a1_phi2; 
    % torque = inertia*a; 

    [t3_phi2, q3_phi2, w3_phi2, torque3_phi2] = gyrostat_discrete_torqueN(dt, t2, t3, inertia, torque_G0, w0, q0); 

    t_phi2 = [t1_phi2; t2_phi2(2:end); t3_phi2(2:end)]; 
    % t_phi2 = [t1_phi2; t1_phi2(end) + t2_phi2(2:end); t1_phi2(end) + t2_phi2(end) + t3_phi2(2:end)]; 
    w_phi2 = [w1_phi2; w2_phi2(2:end ,:); w3_phi2(2:end, :)]; 
    q_phi2 = [q1_phi2; q2_phi2(2:end ,:); q3_phi2(2:end, :)]; 
    torque_phi2 = [torque1_phi2; torque2_phi2(2:end, :); torque3_phi2(2:end, :)]; 

    ypr_phi2 = zeros(length(q_phi2), 3); 
    for i = 1:max(size(q_phi2))
        ypr_phi2(i, :) = SpinCalc('QtoEA321', q_phi2(i, :), eps, 0); 
    end 

    %% Plot phi2
    plot_option = 1; 
    if plot_option == 1
        plot_qwypr(t_phi2, q_phi2, w_phi2, ypr_phi2, 2)
    end 

    time_torque = linspace(t0, t3, length(torque_phi2)); 
    figure
        plot(time_torque, torque_phi2); 
        grid on 
        ylabel('Torque (Nm)') 
        xlabel('Time (sec)') 
        title('Torque (Nm)') 

    %% Determine Phi 3 slew times

    w0 = 0; 
    t0 = 0; 
    wf = 0; 
    [t1, t2, t3] = find_slew_times(aMax, wMax, phi3, w0, wf, t0); 

    % Discretize slew times to 2 decimal places 
    t1 = round(t1, 2); 
    t2 = round(t2, 2); 
    t3 = round(t3, 3); 

    %% Solve for attitue determination - third slew 

    % t0 --> t1 
    tEnd = t1 - t0; 
    w0 = w_phi2(end, :)'; 
    q0 = q_phi2(end, :)'; 
    a = [ 0;  0; aMax];  
    torque = inertia*a; 

    [t1_phi3, y1_phi3] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [0, tEnd], [w0; q0]);

    % t1 --> t2 
    tEnd = t2 - t1; 
    w0 = y1_phi3(end, 1:3)'; 
    q0 = y1_phi3(end, 4:7)'; 
    a = [0; 0; 0]; 
    torque = inertia*a; 

    [t2_phi3, y2_phi3] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [0, tEnd], [w0; q0]); 

    % t2 --> t3 
    tEnd = t3 - t2; 
    w0 = y2_phi3(end, 1:3)'; 
    q0 = y2_phi3(end, 4:7)'; 
    a = [ 0; 0; -aMax]; 
    torque = inertia*a; 

    [t3_phi3, y3_phi3] = ode45(@(t, Z) gyrostat_cont(inertia, torque, Z), [0, tEnd], [w0; q0]); 

    y_phi3 = [y1_phi3; y2_phi3(2:end, :); y3_phi3(2:end, :)]; 
    w_phi3 = y_phi3(:, 1:3); 
    q_phi3 = y_phi3(:, 4:7); 

    ypr_phi3 = zeros(length(q_phi3), 3); 
    for i = 1:max(size(q_phi3))
        ypr_phi3(i, :) = SpinCalc('QtoEA321', q_phi3(i, :), eps, 0); 
    end 

    t_phi3 = [t1_phi3; t1_phi3(end)+t2_phi3(2:end); t1_phi3(end)+t2_phi3(end)+t3_phi3(2:end)]; 

    %% Plot phi3
    plot_option = 0; 
    if plot_option == 1
        plot_qwypr(t_phi3, q_phi3, w_phi3, ypr_phi3, 3)
    end 

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
        
ypr_total = zeros(length(q_total), 3); 
for i = 1:max(size(q_total))
    ypr_total(i, :) = SpinCalc('QtoEA321', q_total(i, :), eps, 0); 
end 
        
%% plot total stuff 

plot_option = 1; 
if plot_option == 1
    plot_qwypr(t_total, q_total, w_total, ypr_total, 'total')
end 
