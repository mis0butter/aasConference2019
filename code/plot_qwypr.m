function plot_qwypr(t, q, w, h, torque, a, ws, phi_num, a_max, w_max, h_max)
% Plot q, w, torque, a (ypr previously) 

fsize = 20; 
phi_num = num2str(phi_num); 

% constraints time 
n = 10; 
t_dt = max(t)/n; 
t_lim = 0:t_dt:max(t); 

% Get y-axis range for plots 
ylimits_q = get_ylimits(q); 
ylimits_w = get_ylimits([-w_max, w_max]); 
% ylimits_ypr = get_ylimits(ypr); 
ylimits_torque = get_ylimits(torque); 
ylimits_a = get_ylimits([-a_max, a_max]); 
ylimits_ws = get_ylimits(ws); 
ylimits_h = get_ylimits([-h_max, h_max]); 

% Angular velocity 
figure('name', strcat(phi_num, '_angular_velocity'))
    plot(t, w(:, 1), 'LineWidth', 1.1)
    hold on; 
    plot(t, w(:, 2), 'k--', 'LineWidth', 1.1)
    plot(t, w(:, 3), 'm-.', 'LineWidth', 1.1)
    plot(t_lim, linspace(w_max, w_max, n+1), 'r:x', t_lim, linspace(-w_max, -w_max,n+1), 'r:x')
    ylim(ylimits_w)
%     ylim([-1.1*w_max, 1.1*w_max]) 
    grid on 
    leg = legend('$\omega_x$', '$\omega_y$', '$\omega_z$', '$\omega_c$', ... 
        'location', 'north', 'orientation', 'horizontal');
        set(leg, 'Interpreter', 'latex', 'FontSize', fsize); 
    ylabel('Angular Velocity (rad/s)', 'Interpreter', 'latex', 'FontSize', fsize) 
    xlabel('Time (sec.)', 'Interpreter', 'latex', 'FontSize', fsize)  
    ax = gca;
    ax.FontSize = fsize - 3;  
%     title(strcat('Angular Velocity: $\dot{\phi}_{', phi_num, '}$'), 'Interpreter', 'latex', 'FontSize', fsize)

% Quaternions 
figure('name', strcat(phi_num, '_quaternions'))
    plot(t, q(:, 1), 'b', 'LineWidth', 1.1)
    hold on; 
    plot(t, q(:, 2), 'r--', 'LineWidth', 1.1)
    plot(t, q(:, 3), 'm-.', 'LineWidth', 1.1)
    plot(t, q(:, 4), 'k:', 'LineWidth', 1.2)
%     plot(t_lim, linspace(1,1,n+1), 'k:x', t_lim, linspace(-1,-1,n+1), 'k:x')
    leg = legend('$q_1$', '$q_2$', '$q_3$', '$q_4$', 'orientation', 'horizontal', 'location', 'north'); 
        set(leg, 'Interpreter', 'latex', 'FontSize', fsize); 
    grid on 
    ylim(ylimits_q)
    ylabel('Quaternions', 'Interpreter', 'latex', 'FontSize', fsize)
    xlabel('Time (sec.)', 'Interpreter', 'latex', 'FontSize', fsize) 
    ax = gca;
    ax.FontSize = fsize - 3;  
%     title(strcat('Quaternion: $q_{', phi_num, '}$'), 'Interpreter', 'latex', 'FontSize', fsize)
    
% Torque 
figure('name', strcat(phi_num, '_torque'))
    plot(t, torque(:, 1), 'LineWidth', 1.1)
    hold on; 
    plot(t, torque(:, 2), 'r--', 'LineWidth', 1.1) 
    plot(t, torque(:, 3), 'm-.', 'LineWidth', 1.1)
%     plot(t_lim, linspace(t_max, t_max, n+1), 'r:x', t_lim, linspace(-t_max, -t_max,n+1), 'r:x')
    leg = legend('$u_x$', '$u_y$', '$u_z$', 'orientation', 'horizontal', 'location', 'north'); 
        set(leg, 'Interpreter', 'latex', 'FontSize', fsize); 
    grid on 
    ylim(ylimits_torque)
%     ylim([-1.1*t_max, 1.1*t_max]) 
    xlabel('Time (sec.)', 'Interpreter', 'latex', 'FontSize', fsize)
    ylabel('Torque (Nm)', 'Interpreter', 'latex', 'FontSize', fsize) 
    ax = gca;
    ax.FontSize = fsize - 3;  
%     title(strcat('Torque: $\tau_{', phi_num, '}$'), 'Interpreter', 'latex', 'FontSize', fsize)
    
% Angular Acceleration 
figure('name', strcat(phi_num, '_angular_acceleration'))
    plot(t(1:end - 1), a(:,1), 'LineWidth', 1.1)
    hold on; 
    plot(t(1:end - 1), a(:, 2), 'k--', 'LineWidth', 1.1)
    plot(t(1:end - 1), a(:, 3), 'm-.', 'LineWidth', 1.1) 
    plot(t_lim, linspace(a_max, a_max, n+1), 'r:x', t_lim, linspace(-a_max, -a_max,n+1), 'r:x')
    leg = legend('$\dot{\omega}_x$', '$\dot{\omega}_y$', '$\dot{\omega}_z$', '$\dot{\omega}_c$', ... 
        'orientation', 'horizontal', 'location', 'north');
        set(leg, 'Interpreter', 'latex', 'FontSize', fsize); 
    grid on 
    ylim(ylimits_a)
%     ylim([-1.1*a_max, 1.1*a_max]) 
    xlabel('Time (sec.)', 'Interpreter', 'latex', 'FontSize', fsize)
    ylabel('Angular Acceleration (rad/s$^2$)', 'Interpreter', 'latex', 'FontSize', fsize) 
    ax = gca;
    ax.FontSize = fsize - 3;  
%     str = strcat('$\omega_{', phi_num, '}$', 'Interpreter', 'latex'); 
%     title(strcat('Angular Acceleration: $\ddot{\phi}_{', phi_num, '}$'), 'Interpreter', 'latex', 'FontSize', fsize)
    
% Angular Momentum 
figure('name', strcat(phi_num, '_angular_momentum'))
    plot(t, h(:,1), 'LineWidth', 1.1)
    hold on; 
    plot(t, h(:, 2), 'k--', 'LineWidth', 1.1)
    plot(t, h(:, 3), 'm-.', 'LineWidth', 1.1) 
    plot(t_lim, linspace(h_max, h_max, n+1), 'r:x', t_lim, linspace(-h_max, -h_max,n+1), 'r:x')
    leg = legend('$h_x$', '$h_y$', '$h_z$', '$h_c$', 'orientation', 'horizontal', 'location', 'north');
        set(leg, 'Interpreter', 'latex', 'FontSize', fsize); 
    grid on 
    ylim(ylimits_h)
%     ylim([-1.1*a_max, 1.1*a_max]) 
    xlabel('Time (sec.)', 'Interpreter', 'latex', 'FontSize', fsize)
    ylabel('Angular Momentum (Nms)', 'Interpreter', 'latex', 'FontSize', fsize) 
    ax = gca;
    ax.FontSize = fsize - 3;  
%     str = strcat('$\omega_{', phi_num, '}$', 'Interpreter', 'latex'); 
%     title(strcat('Angular Acceleration: $\ddot{\phi}_{', phi_num, '}$'), 'Interpreter', 'latex', 'FontSize', fsize)
    
% Wheel Speeds 
figure('name', strcat(phi_num, '_wheel_speeds'))
    plot(t, ws(:,1), 'LineWidth', 1.1)
    hold on; 
    plot(t, ws(:, 2), 'k--', 'LineWidth', 1.1)
    plot(t, ws(:, 3), 'm-.', 'LineWidth', 1.1) 
    plot(t, ws(:, 4), 'm-.', 'LineWidth', 1.1) 
%     plot(t_lim, linspace(h_max, h_max, n+1), 'r:x', t_lim, linspace(-h_max, -h_max,n+1), 'r:x')
    leg = legend('$ws_1$', '$ws_2$', '$ws_3$', '$ws_4$', 'orientation', 'horizontal', 'location', 'north');
        set(leg, 'Interpreter', 'latex', 'FontSize', fsize); 
    grid on 
    ylim(ylimits_ws)
%     ylim([-1.1*a_max, 1.1*a_max]) 
    xlabel('Time (sec.)', 'Interpreter', 'latex', 'FontSize', fsize)
    ylabel('Wheel Speeds (rpm)', 'Interpreter', 'latex', 'FontSize', fsize) 
    ax = gca;
    ax.FontSize = fsize - 3;  
%     str = strcat('$\omega_{', phi_num, '}$', 'Interpreter', 'latex'); 
%     title(strcat('Angular Acceleration: $\ddot{\phi}_{', phi_num, '}$'), 'Interpreter', 'latex', 'FontSize', fsize)
end 
    
% figure()
%     plot(t, ypr(:, 1), 'LineWidth', 1.1)
%     hold on; 
%     plot(t, ypr(:, 2), 'r--', 'LineWidth', 1.1)
%     plot(t, ypr(:, 3), 'g-.', 'LineWidth', 1.1)
%     legend('Yaw', 'Pitch', 'Roll'); 
%     grid on 
%     ylim(ylimits_ypr)
%     xlabel('time (s)') 
%     ylabel('degrees') 
%     title(strcat('Euler Angles: \phi_{', phi_num, '}'))

%% 

function ylimits = get_ylimits(data) 
% Set y limits of axis 

% rng = range(data(:)); 

rng = max(data(:)) - min(data(:)); 

midp = min(data(:)) + rng/2; 
ylimits = [ midp - 1.2*rng/2, midp + 1.6*rng/2 ]; 

end 