function plot_qwypr(t, q, w, ypr, phi_num)
%% Plot q, w, ypr 

phi_num = num2str(phi_num); 

ylimits_q = get_ylimits(q); 
ylimits_w = get_ylimits(w); 
ylimits_ypr = get_ylimits(ypr); 

plot_option = 1; 
if plot_option == 1
figure()
    plot(t, w)
    ylim(ylimits_w)
    legend('w1', 'w2', 'w3'); 
    ylabel('w (rad/s)') 
    xlabel('time (s)') 
    title(strcat('Angular Velocity: Phi=', phi_num))

figure()
    plot(t, q)
    legend('q1', 'q2', 'q3', 'q4'); 
    ylim(ylimits_q)
    ylabel('quats') 
    xlabel('time (s)') 
    title(strcat('Quaternion: Phi=', phi_num)) 
    
figure()
    plot(t, ypr)
    legend('Yaw', 'Pitch', 'Roll'); 
    ylim(ylimits_ypr)
    xlabel('time (s)') 
    ylabel('degrees') 
    title(strcat('Euler Angles: Phi=', phi_num)) 
end 