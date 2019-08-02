function [phi2_P, phi2_P_sum, SP_S] = phi2_curve(S_G0, P1_G0, P2_G0, phi2)
% FINDING PHI2

    SPa_G0 = P2_G0 - S_G0; 
    SPb_G0 = P1_G0 - S_G0; 
    SPe_G0 = cross(SPa_G0, SPb_G0); 

    SP1_G0 = SPb_G0; 
    SP3_G0 = SPe_G0; 
    SP2_G0 = cross(SP1_G0, SP3_G0); 

    SP_DCM_G0 = [SP1_G0'; SP2_G0'; SP3_G0']; 
    G0_DCM_SP = SP_DCM_G0'; 

    dphi = phi2/100;        % radians 

    V_SP(1, :) = [1 0 0];
    V_G0(1, :) = G0_DCM_SP * V_SP(1, :)'; 
    P_G0(1, :) = P1_G0; 
    phi2_P = zeros(99, 1); 
    for i = 1:99 
        V_SP(1 + i, :) = [cos(i*dphi), 0, sin(i*dphi), ]; 
        V_G0(1 + i, :) = G0_DCM_SP * V_SP(1 + i, :)'; 
        P_G0(1 + i, :) = S_G0' + V_G0(1 + i, :); 

        a = P_G0(i, :); 
        b = P_G0(1 + i, :); 
        phi2_P(i) = acos(dot(a, b)/(norm(a)*norm(b))); 
    end 

    phi2_P_sum = sum(phi2_P); 

    plot_option = 1; 
    if plot_option == 1
        figure()
            plot3([0 SP1_G0(1)], [0 SP1_G0(2)], [0 SP1_G0(3)], 'r')
            hold on; grid on 
            plot3([0 SP2_G0(1)], [0 SP2_G0(2)], [0 SP2_G0(3)], 'r')
            plot3([0 SP3_G0(1)], [0 SP3_G0(2)], [0 SP3_G0(3)], 'r')

            plot3([0 Pi_G0(1)], [0 Pi_G0(2)], [0 Pi_G0(3)], 'b'); 
            plot3([0 Pf_G0(1)], [0 Pf_G0(2)], [0 Pf_G0(3)], 'b'); 
            plot3([0 P1_G0(1)], [0 P1_G0(2)], [0 P1_G0(3)], 'b'); 
            plot3([0 P2_G0(1)], [0 P2_G0(2)], [0 P2_G0(3)], 'b'); 
            plot3([0 e_G0(1)], [0 e_G0(2)], [0 e_G0(3)], 'b'); 

            plot3(S_G0 + V_G0(:, 1), S_G0 + V_G0(:, 2), S_G0 + V_G0(:, 3), '-.')
            plot3(P_G0(:, 1), P_G0(:, 2), P_G0(:, 3), '-.')
        %     plot3(GP_G0(:, 1), GP_G0(:, 2), GP_G0(:, 3)) 
    end 
    
end 