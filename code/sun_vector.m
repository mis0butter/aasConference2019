function [alpha, theta_Pi_Sproj, theta_Sproj_Pf, theta_Pi_Pf, S_N, S_PiPf_G0, S_G0] = ... 
    sun_vector(G0_DCM_N, e_G0, Pi_G0, Pf_G0)
% Generate sun vector 

        S_N = [rand*(-1)^round(rand); rand*(-1)^round(rand); rand*(-1)^round(rand)];  
        S_N = S_N/norm(S_N); 
        S_G0 = G0_DCM_N*S_N;               
        alpha = pi/2 - acos(dot(S_G0, e_G0));         % coming out to 0 - check    

        % Sun projection onto slew plane 
        S_PiPf_G0 = cross(e_G0, cross(S_G0, e_G0));
        
%         %% Make alpha = 0  
%         S_G0 = S_PiPf_G0;
%         S_G0 = S_G0/norm(S_G0); 
%         S_N = G0_DCM_N'*S_G0; 
%         alpha = 0;
        
        %% Sun angles with other points 
        
%         S_PiPf_G0 = cross(e_G0, cross(S_G0, e_G0));    % sun projection vector G frame
        S_PiPf_G0 = S_PiPf_G0/norm(S_PiPf_G0);         % sun projection --> unit vector 

        theta_Sproj_Pf = acos(dot(S_PiPf_G0, Pf_G0)); 
        theta_Pi_Pf = acos(dot(Pi_G0, Pf_G0)); 
        theta_Pi_Sproj = acos(dot(Pi_G0, S_PiPf_G0)); 
end 