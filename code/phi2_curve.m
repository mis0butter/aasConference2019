function [phi2_P, phi2_P_sum, SP_S] = phi2_curve(S_G0, P1_G0, P2_G0, phi2)
% FINDING PHI2

S1_G0 = S_G0; 
S2_G0 = P1_G0 - S_G0; 
S3_G0 = -cross(S1_G0, S2_G0); 
S_DCM_G0 = [S1_G0'; S2_G0'; S3_G0']; 
G0_DCM_S = S_DCM_G0'; 

dphi = phi2/100;        % radians 

%%

SP_S(1, :) = [1 0 0];
SP_G0(1, :) = G0_DCM_S * SP_S(1, :)'; 
P_G0(1, :) = P1_G0; 
phi2_P = zeros(99, 1); 
for i = 1:99 
    SP_S(1 + i, :) = [cos(i*dphi), 0, sin(i*dphi), ]; 
    SP_G0(1 + i, :) = G0_DCM_S * SP_S(1 + i, :)'; 
    P_G0(1 + i, :) = S_G0' + SP_G0(1 + i, :); 
    phi2_P(i) = acos(dot(P_G0(i, :), P_G0(1 + i, :))); 
end 

phi2_P_sum = sum(phi2_P); 

%%

plot_option = 1; 
if plot_option == 1
    figure()
        plot3([0 S1_G0(1)], [0 S1_G0(2)], [0 S1_G0(3)])
        hold on; grid on 
        plot3([0 S2_G0(1)], [0 S2_G0(2)], [0 S2_G0(3)])
        plot3([0 S3_G0(1)], [0 S3_G0(2)], [0 S3_G0(3)])

        plot3([0 P1_G0(1)], [0 P1_G0(2)], [0 P1_G0(3)], 'b'); 
        plot3([0 P2_G0(1)], [0 P2_G0(2)], [0 P2_G0(3)], 'b'); 
        plot3(SP_G0(:, 1), SP_G0(:, 2), SP_G0(:, 3), '-.')
        plot3(P_G0(:, 1), P_G0(:, 2), P_G0(:, 3), '-.')
    %     plot3(GP_G0(:, 1), GP_G0(:, 2), GP_G0(:, 3)) 
end 