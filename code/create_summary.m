% Create summary.txt
summary = fopen('summary.txt', 'w'); 

% Write initial, final, sun vector positions 
fprintf(summary, 'Pi in N frame = \t \t \t [%.3f, %.3f, %.3f] \n', Pi_N(1), Pi_N(2), Pi_N(3)); 
fprintf(summary, 'Pf in N frame = \t \t \t [%.3f, %.3f, %.3f] \n', Pf_N(1), Pf_N(2), Pf_N(3)); 
fprintf(summary, 'S in N frame = \t \t \t \t [%.3f, %.3f, %.3f] \n \n', S_N(1), S_N(2), S_N(3)); 

% Write alpha 
fprintf(summary, 'alpha = \t \t \t \t \t %.3f rad, %.3f deg \n \n', alpha, alpha*180/pi); 

% Write v and a constraint s
fprintf(summary, 'aMax = \t \t \t \t \t \t %.3f rad/s^2 \n', aMax); 
fprintf(summary, 'wMax = \t \t \t \t \t \t %.3f rad/s \n \n', wMax); 

% Write phi values 
fprintf(summary, 'phi 1 = \t \t \t \t \t %.3f rad, %.3f deg \n', phi1, phi1*180/pi); 
fprintf(summary, 'phi 2 = \t \t \t \t \t %.3f rad, %.3f deg \n', phi2, phi2*180/pi); 
fprintf(summary, 'phi 3 = \t \t \t \t \t %.3f rad, %.3f deg \n \n', phi2, phi3*180/pi); 

% Print final error 
err_final = acosd(dot(Pf_N, P_phi3_N(end, :)));
fprintf(summary, 'Final error: \t \t \t \t %.2f deg \n', err_final); 

fclose(summary); 
