% post-processing script

%% sun intrusion slew 

% Actual slew profile 
P_G = zeros(length(t_total), 3); 
% convert quaternions --> DCMs --> vectors 
for i = 1:length(t_total)
    G0_DCM_G = quat2DCM(q_total(i, :)); 
    G_DCM_G0 = G0_DCM_G'; 
    P_G(i, :) = G_DCM_G0*Pi_G0; 
end 

% Phi2 slew profile - part 1 (torque) 
P_G1_phi2 = zeros(length(t1_phi2), 3); 
for i = 1:length(t1_phi2)
    G0_DCM_G = quat2DCM(q1_phi2(i, :)); 
    G_DCM_G0 = G0_DCM_G'; 
    P_G1_phi2(i, :) = G_DCM_G0*Pi_G0; 
end 

% Phi2 angle - part 1 (torque) 
for i = 1:length(t1_phi2) - 1
%     a = G0_DCM_G*P_G1_phi2(i, :)'; 
    a = P_G1_phi2(i, :)';           a = a/norm(a); 
    b = S_G0;                       b = b/norm(b); 
    angle1_phi2_SG0(i, 1) = acos(dot(a, b)); 
    b = P_G1_phi2(i + 1, :)';       b = b/norm(b); 
    angle1_phi2_P2G0(i, 1) = acos(dot(a, b)); 
    
end 

% Phi2 slew profile - part 2 (no torque) 
P_G2_phi2 = zeros(length(t2_phi2), 3); 
for i = 1:length(t2_phi2)
    G0_DCM_G = quat2DCM(q2_phi2(i, :)); 
    G_DCM_G0 = G0_DCM_G'; 
    P_G2_phi2(i, :) = G_DCM_G0*Pi_G0; 
end 

% Phi2 angle - part 2 (no torque) 
for i = 1:length(t2_phi2) - 1
%     a = G0_DCM_G*P_G2_phi2(i, :)'; 
    a = P_G2_phi2(i, :)';           a = a/norm(a); 
    b = S_G0;                       b = b/norm(b); 
    angle2_phi2_SG0(i, 1) = acos(dot(a, b)); 
    b = P_G2_phi2(i + 1, :)';       b = b/norm(b); 
    angle2_phi2_P2G0(i, 1) = acos(dot(a, b)); 
end 

% Phi2 slew profile - part 3 (-torque) 
P_G3_phi2 = zeros(length(t3_phi2), 3); 
for i = 1:length(t3_phi2)
    G0_DCM_G = quat2DCM(q3_phi2(i, :)); 
    G_DCM_G0 = G0_DCM_G'; 
    P_G3_phi2(i, :) = G_DCM_G0*Pi_G0; 
end 

% Phi2 angle - part 3 (-torque) 
for i = 1:length(t3_phi2) - 1
%     a = G0_DCM_G*P_G3_phi2(i, :)'; 
    a = P_G3_phi2(i, :)';           a = a/norm(a); 
    b = S_G0;                       b = b/norm(b); 
    angle3_phi2_SG0(i, 1) = acos(dot(a, b)); 
    b = P_G3_phi2(i + 1, :)';       b = b/norm(b); 
    angle3_phi2_P2G0(i, 1) = acos(dot(a, b)); 
end 



%% NO sun intrusion slew - nominal 

P_Gnom = zeros(length(t_phiNom), 3); 
% convert quaternions --> DCMs --> vectors 
for i = 1:length(t_phiNom)
    G0_DCM_G = quat2DCM(q_phiNom(i, :)); 
    G_DCM_G0 = G0_DCM_G'; 
    P_Gnom(i, :) = G_DCM_G0*Pi_G0; 
end 

% Phi3 angle - 

%%

figure()

    % Actual Slew Profile 
    plot3(P_G(:, 1), P_G(:, 2), P_G(:,3), 'r', 'LineWidth', 1.1); 
    
    grid on; hold on
    
%     % Phi2 Slew Profile 
%     plot3(P_G1_phi2(:, 1), P_G1_phi2(:, 2), P_G1_phi2(:, 3), 'r'); 
%     plot3(P_G2_phi2(:, 1), P_G2_phi2(:, 2), P_G2_phi2(:, 3), 'g');
%     plot3(P_G3_phi2(:, 1), P_G3_phi2(:, 2), P_G3_phi2(:, 3), 'r');
    
    % Nominal (no sun intrusion) profile 
%     plot3(P_Gnom(:, 1), P_Gnom(:, 2), P_Gnom(:,3), '-.'); 

    % Phi1 and 3 profiles - taken from nominal profile 
%     plot3(phi1_P(:, 1), phi1_P(:, 2), phi1_P(:, 3), '-.'); 
%     plot3(P_phi3_G0(:, 1), P_phi3_G0(:, 2), P_phi3_G0(:, 3), '-.'); 
%     plot3(P_phi1_G0(:, 1), P_phi1_G0(:, 2), P_phi1_G0(:, 3), '-.'); 

    % desired phi2 slew profile 
%     plot3(P_phi2_G0(:, 1), P_phi2_G0(:, 2), P_phi2_G0(:, 3), '-.')
    
    plot3([0 Pi_G0(1)], [0 Pi_G0(2)], [0 Pi_G0(3)], 'b:', 'LineWidth', 1.1); 
    plot3([0 Pf_G0(1)], [0 Pf_G0(2)], [0 Pf_G0(3)], 'b:', 'LineWidth', 1.1); 
    plot3([0 P1_G0(1)], [0 P1_G0(2)], [0 P1_G0(3)], 'b:', 'LineWidth', 1.1); 
    plot3([0 P2_G0(1)], [0 P2_G0(2)], [0 P2_G0(3)], 'b:', 'LineWidth', 1.1); 
    plot3([0 e_G0(1)*0.5], [0 e_G0(2)*0.5], [0 e_G0(3)*0.5], 'g-.'); 
     
%     plot3([0 S1_G0(1)], [0 S1_G0(2)], [0 S1_G0(3)])
%     plot3([0 S2_G0(1)], [0 S2_G0(2)], [0 S2_G0(3)])
%     plot3([0 S3_G0(1)], [0 S3_G0(2)], [0 S3_G0(3)])

%     plot3(P_S(:, 1), P_S(:, 2), P_S(:, 3))
    
    plot(theta,sin(theta)./theta,'LineWidth',3) 
    
    % Sun vector 
    plot3(S_G0(1), S_G0(2), S_G0(3), 'p'); 
    % Sun projection 
    plot3([0 S_PiPf_G0(1)], [0 S_PiPf_G0(2)], [0 S_PiPf_G0(3)], 'g-.'); 
    
    % Phi2 - P3 (P1,P2 proj onto Sun vector) lines 
    plot3([0 P3_G0(1)], [0 P3_G0(2)], [0 P3_G0(3)], 'b:', 'LineWidth', 1.1); 
    plot3([P3_G0(1) P1_G0(1)], [P3_G0(2) P1_G0(2)], [P3_G0(3) P1_G0(3)], 'r:'); 
    plot3([P3_G0(1) P2_G0(1)], [P3_G0(2) P2_G0(2)], [P3_G0(3) P2_G0(3)], 'r:'); 
    
    % Initial, final, P2, P2 vectors 
    text(Pi_G0(1), Pi_G0(2), Pi_G0(3), sprintf(' P_i')) 
    text(Pf_G0(1), Pf_G0(2), Pf_G0(3), sprintf(' P_f')) 
    text(e_G0(1)*0.5, e_G0(2)*0.5, e_G0(3)*0.5, sprintf(' e')) 
    text(P1_G0(1), P1_G0(2), P1_G0(3), sprintf(' P_1')) 
    text(P2_G0(1), P2_G0(2), P2_G0(3), sprintf(' P_2')) 
    text(S_G0(1), S_G0(2), S_G0(3), sprintf(' sun')) 
%     text(S_PiPf_G0(1), S_PiPf_G0(2), S_PiPf_G0(3), sprintf(' sun_{proj}')) 
    
    xlabel('G0_x')
    ylabel('G0_y') 
    zlabel('G0_z') 
    title('\phi_1, \phi_2, and \phi_3 Slews') 
    
    
    