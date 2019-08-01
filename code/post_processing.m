% post-processing script

%% sun intrusion slew 

P_G = zeros(length(t_total), 3); 
% convert quaternions --> DCMs --> vectors 
for i = 1:length(t_total)
    G0_DCM_G = quat2DCM(q_total(i, :)); 
    G_DCM_G0 = G0_DCM_G'; 
    P_G(i, :) = G_DCM_G0*Pi_G0; 
end 

%% NO sun intrusion slew - nominal 

P_Gnom = zeros(length(t_phiNom), 3); 
% convert quaternions --> DCMs --> vectors 
for i = 1:length(t_phiNom)
    G0_DCM_G = quat2DCM(q_phiNom(i, :)); 
    G_DCM_G0 = G0_DCM_G'; 
    P_Gnom(i, :) = G_DCM_G0*Pi_G0; 
end 

%%

figure()
    plot3(P_G(:, 1), P_G(:, 2), P_G(:,3)); 
    grid on; hold on
%     plot3(P_Gnom(:, 1), P_Gnom(:, 2), P_Gnom(:,3)); 
    
    plot3([0 Pi_G0(1)], [0 Pi_G0(2)], [0 Pi_G0(3)], 'b'); 
    plot3([0 Pf_G0(1)], [0 Pf_G0(2)], [0 Pf_G0(3)], 'b'); 
    plot3([0 P1_G0(1)], [0 P1_G0(2)], [0 P1_G0(3)], 'b'); 
    plot3([0 P2_G0(1)], [0 P2_G0(2)], [0 P2_G0(3)], 'b'); 
    plot3([0 e_G0(1)], [0 e_G0(2)], [0 e_G0(3)], 'b'); 
    
    plot(theta,sin(theta)./theta,'LineWidth',3) 
    
    plot3([0 S_G0(1)], [0 S_G0(2)], [0 S_G0(3)], 'r'); 
    plot3([0 S_PiPf_G0(1)], [0 S_PiPf_G0(2)], [0 S_PiPf_G0(3)], 'r'); 
    plot3([0 P3_G0(1)], [0 P3_G0(2)], [0 P3_G0(3)], 'g'); 
    
    plot3([P3_G0(1) P1_G0(1)], [P3_G0(2) P1_G0(2)], [P3_G0(3) P1_G0(3)], 'g'); 
    plot3([P3_G0(1) P2_G0(1)], [P3_G0(2) P2_G0(2)], [P3_G0(3) P2_G0(3)], 'g'); 
    
    text(P3_G0(1), P3_G0(2), P3_G0(3), ... 
        sprintf('    phi2 = %0.2f deg', phi2_P3*180/pi))
    text(Pi_G0(1), Pi_G0(2), Pi_G0(3), sprintf(' Pi')) 
    text(Pf_G0(1), Pf_G0(2), Pf_G0(3), sprintf(' Pf')) 
    text(e_G0(1), e_G0(2), e_G0(3), sprintf(' e')) 
    text(P1_G0(1), P1_G0(2), P1_G0(3), sprintf(' P1')) 
    text(P2_G0(1), P2_G0(2), P2_G0(3), sprintf(' P2')) 
    text(S_G0(1), S_G0(2), S_G0(3), sprintf(' sun')) 
    text(S_PiPf_G0(1), S_PiPf_G0(2), S_PiPf_G0(3), sprintf(' sun proj')) 
    
    xlabel('G0_x')
    ylabel('G0_y') 
    zlabel('G0_z') 
    
    
    