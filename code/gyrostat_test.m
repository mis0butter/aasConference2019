% gyrostat test functions 

w0 = [-0.1; -0.2; 0.3]; 
q0 = [cosd(45); 0; 0; cosd(45)]; 
a = [-1; 0.5; -0.3]; 
inertia = [ 100     0       0; 
            0       100     0; 
            0       0       100]; 
torque0 = inertia*a; 
tEnd = 10; 

%%

% Continuous applies constant torque in the BODY frame 
[t_c, y_c] = ode45(@(t,Z) gyrostat_cont(inertia, torque0, Z), [0, tEnd], [w0; q0]);
w_c = y_c(:, 1:3); 
q_c = y_c(:, 4:7); 

dt = 1/100; 

% Following discrete solvers apply constant torque in the INERTIAL frame 
[t_d, q_d, w_d, torque_d] = gyrostat_discrete(dt, 0, tEnd, inertia, torque0, w0, q0);               % torque in --> body frame 
[t_cd, q_cd, w_cd, torque_cd] = gyrostat_CD(dt, 0, tEnd, inertia, torque0, w0, q0);                 % torque in --> body frame 
[t_cdN, q_cdN, w_cdN, torque_cdN] = gyrostat_CD_torqueN(dt, 0, tEnd, inertia, torque0, w0, q0);     % torque in --> inertial frame 
[t_dN, q_dN, w_dN, torque_dN] = gyrostat_discrete_torqueN(dt, 0, tEnd, inertia, torque0, w0, q0);   % torque in --> inertial frame 
    
%% Plots 

figure
for i = 1:4
    subplot(4,1,i) 
        plot(t_c, q_c(:,i), t_d, q_d(:,i), t_cd, q_cd(:, i), '--', ... 
            t_dN, q_dN(:, i), '-.', t_cdN, q_cdN(:, i), ':')
        grid on
        ylabel(strcat('q', num2str(i)))
        if i == 1
            title('Continuous and Discrete Quaternions') 
            legend('Continuous', 'Discrete', 'Continuous-Discrete', ... 
                'Discrete torqueN', 'Continuous-Discrete torqueN') 
        end 
end 
xlabel('Time (sec)') 

figure
for i = 1:3
    subplot(3,1,i)
        plot(t_c, w_c(:, i), t_d, w_d(:, i), t_cd, w_cd(:, i), '--', ... 
            t_dN, w_dN(:, i), '-.', t_cdN, w_cdN(:, i), ':')
        grid on 
        ylabel(strcat('w', num2str(i)))
        if i == 1
            title('Continuous and Discrete Angular Velocity') 
            legend('Continuous', 'Discrete', 'Continuous-Discrete', ... 
                'Discrete torqueN', 'Continuous-Discrete torqueN') 
        end 
end 
xlabel('Time (sec)')