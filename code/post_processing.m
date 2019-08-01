% post-processing script 

P_G = zeros(length(t_total), 3); 
% convert quaternions --> DCMs --> vectors 
for i = 1:length(t_total)
    G0_DCM_G = quat2DCM(q_total(i, :)); 
    G_DCM_G0 = G0_DCM_G'; 
    P_G(i, :) = G_DCM_G0*Pi_G0; 
end 

%%

figure()
    plot3(P_G(:, 1), P_G(:, 2), P_G(:,3)); 
    grid on; hold on
    xlabel('G0_x')
    ylabel('G0_y') 
    zlabel('G0_z') 
    
    
    