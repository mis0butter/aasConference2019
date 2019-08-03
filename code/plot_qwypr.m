function plot_qwypr(t, q, w, torque, ypr, a, phi_num)
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
    plot(t, w)
    ylim(ylimits_w)
    grid on 
    legend('w1', 'w2', 'w3'); 
    ylabel('w (rad/s)') 
    xlabel('time (s)') 
    title(strcat('Angular Velocity: Phi=', phi_num))

figure()
    plot(t, q)
    legend('q1', 'q2', 'q3', 'q4'); 
    grid on 
    ylim(ylimits_q)
    ylabel('quats') 
    xlabel('time (s)') 
    title(strcat('Quaternion: Phi=', phi_num)) 
    
figure()
    plot(t, ypr)
    legend('Yaw', 'Pitch', 'Roll'); 
    grid on 
    ylim(ylimits_ypr)
    xlabel('time (s)') 
    ylabel('degrees') 
    title(strcat('Euler Angles: Phi=', phi_num)) 
    
figure()
    plot(t, torque)
    legend('Gx', 'Gy', 'Gz'); 
    grid on 
    ylim(ylimits_torque)
    xlabel('time (s)') 
    ylabel('Nm') 
    title(strcat('Torque: Phi=', phi_num)) 
    
figure()
    plot(t(1:end - 1), a)
    legend('Gx', 'Gy', 'Gz'); 
    grid on 
    ylim(ylimits_a)
    xlabel('time (s)') 
    ylabel('rad/s^2') 
    title(strcat('Angular Acceleration: Phi=', phi_num)) 
end 