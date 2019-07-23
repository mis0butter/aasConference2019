% gyrostat test functions 

w0 = [-1; -2; 3]; 
q0 = [0; 0; 0; 1]; 
a = [-1; 0.5; -0.3]; 
inertia = [ 100     0       0; 
            0       100     0; 
            0       0       100]; 
torque0 = inertia*a; 
tEnd = 5; 

%%

% Continuous may NOT be correct - direction of torque changes with each
% step. 
[t_c, y_c] = ode45(@(t,Z) gyrostat_cont(inertia, torque0, Z), [0, tEnd], [w0; q0]);
w_c = y_c(:, 1:3); 
q_c = y_c(:, 4:7); 

dt = 1/100; 

[t_d, q_d, w_d, torque_d] = gyrostat_discrete(dt, 0, tEnd, inertia, torque0, w0, q0); 

[t_cd, q_cd, w_cd, torque_cd] = gyrostat_CD(dt, 0, tEnd, inertia, torque0, w0, q0); 

[t_dN, q_dN, w_dN, torque_dN] = gyrostat_discrete_torqueN(dt, 0, tEnd, inertia, torque0, w0, q0); 
    
%% Plots 

figure
for i = 1:4
    subplot(4,1,i) 
        plot(t_c, q_c(:,i), t_d, q_d(:,i), '-.', t_cd, q_cd(:, i), '--', t_dN, q_dN(:, i), ':')
        grid on
        ylabel(strcat('q', num2str(i)))
        if i == 1
            title('Continuous and Discrete Quaternions') 
            legend('Continuous', 'Discrete', 'Continuous-Discrete', 'Discrete torqueN') 
        end 
end 
xlabel('Time (sec)') 

figure
for i = 1:3
    subplot(3,1,i)
        plot(t_c, w_c(:, i), t_d, w_d(:, i), '-.', t_cd, w_cd(:, i), '--', t_dN, w_dN(:, i), ':')
        grid on 
        ylabel(strcat('w', num2str(i)))
        if i == 1
            title('Continuous and Discrete Angular Velocity') 
            legend('Continuous', 'Discrete', 'Continuous-Discrete', 'Discrete torqueN') 
        end 
end 
xlabel('Time (sec)')