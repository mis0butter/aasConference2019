% SP_plane 

%%%

% %%
% 
% % SP plane is not actually from sun vector to P --> it is now P3 to P
% % plane, but we are keeping the nomenclature. Sorry. 
% 
% % FINDING PHI2
% 
%     SPa_G0 = P2_G0 - P3_G0; 
%     SPb_G0 = P1_G0 - P3_G0; 
%     SPe_G0 = cross(SPa_G0, SPb_G0); 
% 
%     SP1_G0 = SPb_G0; 
%         SP1_G0 = SP1_G0 / norm(SP1_G0); 
%     SP3_G0 = SPe_G0; 
%         SP3_G0 = SP3_G0 / norm(SP3_G0); 
%     SP2_G0 = cross(SP1_G0, SP3_G0); 
%         
%     SP_DCM_G0 = [SP1_G0'; SP2_G0'; SP3_G0']; 
%     G0_DCM_SP = SP_DCM_G0'; 
%     
%     P3_SP = SP_DCM_G0*P3_G0; 
% 
%     dphi = 0.001;        % radians 
% 
%     v_length = norm(P1_G0 - P3_G0); 
%     V_SP(1, :) = [1 0 0]*v_length;
%     P_SP(1, :) = SP_DCM_G0*P1_G0; 
%     V_G0(1, :) = G0_DCM_SP * V_SP(1, :)'; 
%     P_phi2_G0(1, :) = P1_G0; 
%     
%     i = 1; 
%     P_phi2_G0(i, :) = P1_G0'; 
%     while dot(P_phi2_G0(i, :), P2_G0') < 0.9999
%         
%         i = i + 1; 
%         V_SP(i, :) = [cos(i*dphi), sin(i*dphi), 0 ] * v_length; 
%         P_SP(i, :) = P3_SP' + V_SP(i, :); 
%         
%         V_G0(i, :) = G0_DCM_SP * V_SP(i, :)'; 
%         P_phi2_G0(i, :) = P3_G0' + V_G0(i, :); 
%         P_phi2_G0(i, :) = P_phi2_G0(i, :) / norm(P_phi2_G0(i, :)); 
% 
%         a = P_phi2_G0(i, :); 
%         b = P_phi2_G0(i - 1, :); 
%         phi2_P(i) = acos(dot(a, b)/(norm(a)*norm(b)));  
%         
%     end 
% 
%     phi2_P_sum = sum(phi2_P); 
    