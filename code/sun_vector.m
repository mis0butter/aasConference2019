function [alpha, theta_Pi_Sproj, theta_Sproj_Pf, theta_Pi_Pf, S_N, S_PiPf_N] = ... 
    sun_vector(e_N, Pi_N, Pf_N)
% Generate sun vector 

    S_N = [rand*(-1)^round(rand); rand*(-1)^round(rand); rand*(-1)^round(rand)];  
    S_N = S_N/norm(S_N);              
    alpha = pi/2 - acos(dot(S_N, e_N));         % coming out to 0 - check    

    % Sun projection onto slew plane 
    S_PiPf_N = cross(e_N, cross(S_N, e_N));

%         %% Make alpha = 0  
%         S_G0 = S_PiPf_G0;
%         S_G0 = S_G0/norm(S_G0); 
%         S_N = G0_DCM_N'*S_G0; 
%         alpha = 0;

    %% Sun angles with other points 

%     S_PiPf_G0 = cross(e_G0, cross(S_G0, e_G0));    % sun projection vector G frame
    S_PiPf_N = S_PiPf_N/norm(S_PiPf_N);         % sun projection --> unit vector 

    theta_Sproj_Pf = acos(dot(S_PiPf_N, Pf_N)); 
    theta_Pi_Pf = acos(dot(Pi_N, Pf_N)); 
    theta_Pi_Sproj = acos(dot(Pi_N, S_PiPf_N)); 
end 