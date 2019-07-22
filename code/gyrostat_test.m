% gyrostat test functions 

w_in = [0.5; 0; 0]; 
q_in = [0; 0; 0; 1]; 
torque_x = [1; 2; 0.1]; 
inertia = [ 408     0       0; 
            0       427     0; 
            0       0       305]; 
tEnd = 20; 

% Continuous may not be correct - direction of torque changes with each
% step. 
[t_c, y_c] = ode45(@(t,Z) gyrostat_cont(inertia, torque_x, Z), [0, tEnd], [w_in; q_in]);
w_c = y_c(:, 1:3); 
q_c = y_c(:, 4:7); 

dt = 1/100; 

[t_d, q_d, w_d] = gyrostat_discrete(dt, 0, tEnd, inertia, torque_x, w_in, q_in); 

%% Plots 

figure
for i = 1:4
    subplot(4,1,i) 
        plot(t_c, q_c(:,i), t_d, q_d(:,i))
        grid on
        ylabel(strcat('q', num2str(i)))
        if i == 1
            title('Continuous and Discrete Quaternions') 
            legend('Continuous', 'Discrete') 
        end 
end 
xlabel('Time (sec)') 

figure
for i = 1:3
    subplot(3,1,i)
        plot(t_c, w_c(:, i), t_d, w_d(:, i))
        grid on 
        ylabel(strcat('w', num2str(i)))
        if i == 1
            title('Continuous and Discrete Angular Velocity') 
            legend('Continuous', 'Discrete') 
        end 
end 
xlabel('Time (sec)')