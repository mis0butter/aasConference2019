function plot_qwypr(t, q, w, torque, a, ypr, phi_num)
% Plot q, w, torque, ypr 

phi_num = num2str(phi_num); 

ylimits_q = get_ylimits(q); 
ylimits_w = get_ylimits(w); 
ylimits_ypr = get_ylimits(ypr); 
ylimits_torque = get_ylimits(torque); 
ylimits_a = get_ylimits(a); 

plot_option = 1; 
if plot_option == 1
figure()
    plot(t, w(:, 1), 'LineWidth', 1.1)
    hold on; 
    plot(t, w(:, 2), 'r--', 'LineWidth', 1.1)
    plot(t, w(:, 3), 'g-.', 'LineWidth', 1.1)
    ylim(ylimits_w)
    grid on 
    legend('\omega_1', '\omega_2', '\omega_3'); 
    ylabel('\omega (rad/s)') 
    xlabel('time (s)') 
    title(strcat('Angular Velocity: \omega_{', phi_num, '}'))

figure()
    plot(t, q(:, 1), 'b', 'LineWidth', 1.1)
    hold on; 
    plot(t, q(:, 2), 'r--', 'LineWidth', 1.1)
    plot(t, q(:, 3), 'g-.', 'LineWidth', 1.1)
    plot(t, q(:, 4), 'k:', 'LineWidth', 1.2)
    legend('q_1', 'q_2', 'q_3', 'q_4'); 
    grid on 
    ylim(ylimits_q)
    ylabel('quats') 
    xlabel('time (s)') 
    title(strcat('Quaternion: q_{', phi_num, '}'))
    
figure()
    plot(t, ypr(:, 1), 'LineWidth', 1.1)
    hold on; 
    plot(t, ypr(:, 2), 'r--', 'LineWidth', 1.1)
    plot(t, ypr(:, 3), 'g-.', 'LineWidth', 1.1)
    legend('Yaw', 'Pitch', 'Roll'); 
    grid on 
    ylim(ylimits_ypr)
    xlabel('time (s)') 
    ylabel('degrees') 
    title(strcat('Euler Angles: \phi_{', phi_num, '}'))
    
figure()
    plot(t, torque(:, 1), 'LineWidth', 1.1)
    hold on; 
    plot(t, torque(:, 2), 'r--', 'LineWidth', 1.1) 
    plot(t, torque(:, 3), 'g-.', 'LineWidth', 1.1)
    legend('Gx', 'Gy', 'Gz'); 
    grid on 
    ylim(ylimits_torque)
    xlabel('time (s)') 
    ylabel('Nm') 
    title(strcat('Torque_{', phi_num, '}'))
    
figure()
    plot(t(1:end - 1), a(:,1), 'LineWidth', 1.1)
    hold on; 
    plot(t(1:end - 1), a(:, 2), 'r--', 'LineWidth', 1.1)
    plot(t(1:end - 1), a(:, 3), 'g-.', 'LineWidth', 1.1) 
    legend('\alpha_x', '\alpha_y', '\alpha_z'); 
    grid on 
    ylim(ylimits_a)
    xlabel('time (s)') 
    ylabel('rad/s^2') 
    title(strcat('Angular Acceleration: \alpha_{', phi_num, '}'))
end 