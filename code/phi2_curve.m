% FINDING PHI2

S2_G0 = S_G0; 
S1_G0 = P1_G0 - S_G0; 
S3_G0 = cross(S1_G0, S2_G0); 
S_DCM_G0 = [S1_G0'; S2_G0'; S3_G0']; 
G0_DCM_S = S_DCM_G0'; 

dphi = phi2/100; 

%%

P_S(1, :) = (S_DCM_G0*P1_G0)'; 
P_G0(1, :) = P1_G0; 
phi2 = zeros(99, 1); 
for i = 1:99 
    P_S(1 + i, :) = [cos(i*dphi), 0, sin(i*dphi)]; 
    P_G0(1 + i, :) = G0_DCM_S*P_S(i, :)'; 
    phi2(i) = acos(dot(P_G0(i, :), P_G0(1 + i, :))); 
end 

phi2 = sum(phi2); 

%%

figure()
    plot3([0 S1_G0(1)], [0 S1_G0(2)], [0 S1_G0(3)])
    hold on; grid on 
    plot3([0 S2_G0(1)], [0 S2_G0(2)], [0 S2_G0(3)])
    plot3([0 S3_G0(1)], [0 S3_G0(2)], [0 S3_G0(3)])
    plot3(P_G0(:, 1), P_G0(:, 2), P_G0(:, 3))
%     plot3(GP_G0(:, 1), GP_G0(:, 2), GP_G0(:, 3)) 