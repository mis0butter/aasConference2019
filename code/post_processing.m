% post-processing script
% Junette Hsin 
    
%% Plot q, w, torque, a

% phi1
plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_phi1, q_phi1, w_phi1, torque_phi1, a_phi1, 'phi1', aMax, wMax)
end 

% phi2
plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_phi2, q_phi2, w_phi2, torque_phi2, a_phi2, 'phi2', aMax, wMax)
end 

% phi3
plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_phi3, q_phi3, w_phi3, torque_phi3, a_phi3, 'phi3', aMax, wMax)
end 

% phi total 
plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_total, q_total, w_total, torque_total, a_total, 'total', aMax, wMax)
end 

% phi nominal 
plot_option = 0; 
if plot_option == 1
    plot_qwypr(t_phiNom, q_phiNom, w_phiNom, torque_phiNom, a_phiNom, 'phiNom', aMax, wMax)
end 

%% quats --> DCMs --> vectors 

% Actual slew profile 
P_G = zeros(length(t_total), 3); 
P_N = zeros(length(t_total), 3); 
% convert quaternions --> DCMs --> vectors 
for i = 1:length(t_total)
    G0_DCM_G = quat2DCM(q_total(i, :)); 
    G_DCM_G0 = G0_DCM_G'; 
    P_G(i, :) = G_DCM_G0*Pi_G0; 
    P_N(i, :) = N_DCM_G0*P_G(i, :)'; 
end 

% phi1 slew 
P_phi1_N = zeros(length(q_phi1), 3); 
for i = 1:length(t_phi1) 
    G0_DCM_G = quat2DCM(q_phi1(i, :)); 
    P_phi1_N(i, :) = N_DCM_G0*G0_DCM_G'*Pi_G0; 
end 

% phi12 slew 
P_phi2_N = zeros(length(q_phi2), 3); 
for i = 1:length(t_phi2) 
    G0_DCM_G = quat2DCM(q_phi2(i, :)); 
    P_phi2_N(i, :) = N_DCM_G0*G0_DCM_G'*Pi_G0; 
end 

% phi3 slew 
P_phi3_N = zeros(length(q_phi3), 3); 
for i = 1:length(t_phi3) 
    G0_DCM_G = quat2DCM(q_phi3(i, :)); 
    P_phi3_N(i, :) = N_DCM_G0*G0_DCM_G'*Pi_G0; 
end 

%% NO sun intrusion slew - nominal 

P_nom_G = zeros(length(t_phiNom), 3); 
P_nom_N = zeros(length(t_phiNom), 3); 
% convert quaternions --> DCMs --> vectors 
for i = 1:length(t_phiNom)
    G0_DCM_G = quat2DCM(q_phiNom(i, :)); 
    G_DCM_G0 = G0_DCM_G'; 
    P_nom_G(i, :) = G_DCM_G0*Pi_G0; 
    P_nom_N(i, :) = N_DCM_G0*P_nom_G(i, :)'; 
end 

%% PLOT!

fsize = 14; 

plot_option = 0; 
if plot_option == 1
    
    figure
    
    % Transparent sphere 
    [x, y, z] = sphere(100);
    h = surfl(x, y, z); 
    set(h, 'FaceAlpha', 0.1)
    shading interp
    
    grid on; hold on
    
    % Plot phi1, phi2, phi3 slew 
%     plot3(P_N(:, 1), P_N(:, 2), P_N(:,3), 'r', 'LineWidth', 2); 
    plot3(P_phi1_N(:, 1), P_phi1_N(:, 2), P_phi1_N(:,3), 'm', 'LineWidth', 2); 
    plot3(P_phi2_N(:, 1), P_phi2_N(:, 2), P_phi2_N(:,3), 'c', 'LineWidth', 2); 
    plot3(P_phi3_N(:, 1), P_phi3_N(:, 2), P_phi3_N(:,3), 'y', 'LineWidth', 2); 
    
%     % Plot nominal slew 
%     plot3(P_nom_N(:, 1), P_nom_N(:, 2), P_nom_N(:, 3), 'b', 'LineWidth', 2); 
    
    % Plot Pi, Pf, P1, P2 points 
    plot3([0 Pi_N(1)], [0 Pi_N(2)], [0 Pi_N(3)], 'b:', 'LineWidth', 1.1); 
        plot3(Pi_N(1), Pi_N(2), Pi_N(3), 'ko', 'LineWidth', 1)
    plot3([0 Pf_N(1)], [0 Pf_N(2)], [0 Pf_N(3)], 'b:', 'LineWidth', 1.1); 
        plot3(Pf_N(1), Pf_N(2), Pf_N(3), 'ko', 'LineWidth', 1)
    plot3([0 P1_N(1)], [0 P1_N(2)], [0 P1_N(3)], 'b:', 'LineWidth', 1.1); 
        plot3(P1_N(1), P1_N(2), P1_N(3), 'ko', 'LineWidth', 1)
    plot3([0 P2_N(1)], [0 P2_N(2)], [0 P2_N(3)], 'b:', 'LineWidth', 1.1); 
        plot3(P2_N(1), P2_N(2), P2_N(3), 'ko', 'LineWidth', 1)
    plot3([0 e_N(1)], [0 e_N(2)], [0 e_N(3)], 'g-.', 'LineWidth', 1.1); 
    
%     plot(theta,sin(theta)./theta,'LineWidth',3) 
    
    % Sun vector 
    plot3(S_N(1), S_N(2), S_N(3), 'rp', 'LineWidth', 3); 
%     % Sun projection 
%     plot3([0 S_PiPf_N(1)], [0 S_PiPf_N(2)], [0 S_PiPf_N(3)], 'g-.'); 
    
    % Phi2 - P3 (P1,P2 proj onto Sun vector) lines 
    plot3([0 P3_N(1)], [0 P3_N(2)], [0 P3_N(3)], 'b:', 'LineWidth', 1.1); 
    plot3([P3_N(1) P1_N(1)], [P3_N(2) P1_N(2)], [P3_N(3) P1_N(3)], 'r:', 'LineWidth', 1.1); 
    plot3([P3_N(1) P2_N(1)], [P3_N(2) P2_N(2)], [P3_N(3) P2_N(3)], 'r:', 'LineWidth', 1.1); 
    
    % Initial, final, P2, P2 vectors 
    t(2) = text(Pi_N(1), Pi_N(2), Pi_N(3), strcat('\phantom{   }', sprintf('$P_i$')), 'Interpreter', 'latex', 'FontSize', fsize); 
    t(2) = text(Pf_N(1), Pf_N(2), Pf_N(3), strcat('\phantom{   }', sprintf('$P_f$')), 'Interpreter', 'latex', 'FontSize', fsize); 
    t(3) = text(e_N(1), e_N(2), e_N(3), strcat('\phantom{   }', sprintf('e')), 'Interpreter', 'latex', 'FontSize', fsize); 
    t(4) = text(P1_N(1), P1_N(2), P1_N(3), strcat('\phantom{   }', sprintf('$P_1$')), 'Interpreter', 'latex', 'FontSize', fsize); 
    t(5) = text(P2_N(1), P2_N(2), P2_N(3), strcat('\phantom{   }', sprintf('$P_2$')), 'Interpreter', 'latex', 'FontSize', fsize); 
    t(6) = text(S_N(1), S_N(2), S_N(3), strcat('\phantom{   }', sprintf('sun')), 'Interpreter', 'latex', 'FontSize', fsize); 
%     text(S_PiPf_N(1), S_PiPf_N(2), S_PiPf_N(3), sprintf(' sun_{proj}')) 
    
    xlabel('$N_x$','interpreter','latex', 'FontSize', fsize)
    ylabel('$N_y$','interpreter','latex', 'FontSize', fsize)
    zlabel('$N_z$','interpreter','latex', 'FontSize', fsize)
    ax = gca;
    ax.FontSize = fsize - 2;  
%     title('$\phi_1$, $\phi_2$, and $\phi_3$ Slews', 'Interpreter', 'latex', 'FontSize', fsize)

end 

%% how far off? 

phi_total = 180/pi*(phi1 + phi2 + phi3); 
disp(strcat(sprintf('phi 1+2+3: %.2f', phi_total), ' deg'))

err_final = acosd(dot(Pf_N, P_phi3_N(end, :)));
disp(strcat(sprintf('Final error: %.2f', err_final), ' deg'))
