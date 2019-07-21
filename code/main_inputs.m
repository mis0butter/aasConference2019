

%% Inputs 

% Non-changing parameters; should be function inputs 
inertia = [ 408     0       0; 
            0       427     0; 
            0       0       305]; 
Pi_G = [0; 1; 0];                           % Pi = unit vector of initial point in the G frame 
Pf_G = [1; 0; 0];                           % Pf = unit vector of the final point in the G frame 
S_N = [cosd(45); cosd(45); cosd(45)];       % S = unit vector of sun vector in the N frame 
S_N = S_N/norm(S_N);                        % normalizing sun vector 
% G_DCM_N = angle2dcm(1, 1, 1);         % DCM from the N to the G frame 
G0_DCM_N = eye(3);                           % DCM from the N to the G frame 
N_DCM_G0 = G0_DCM_N';                         % G to N frame - initial!!! G frame will change throughout sim 
S_G = G0_DCM_N*S_N;                          % Sun vector in G frame 
ep = pi/12;                                 % payload half-cone angle. pi/12 rad = 15 deg  
aMax = 1;                                  % Maximum acceleration, rad/s^2
wMax = 1;                                  % Maximum angular velocity, rad/s

%% Calculate slew angles 

% Calculate normal vector of slew plane 
e_G = cross(Pf_G, Pi_G) / norm(cross(Pf_G, Pi_G));  % eigenaxis of PiPf plane, in G frame 
Pperp_G = cross(Pi_G, e_G);              % perpendicular vector to e and Pi, slew plane, G frame 
P_DCM_G = [Pi_G'; Pperp_G'; e_G'];     % DCM from G frame to P (slew plane) frame 
G_DCM_P = P_DCM_G';                         % P to G frame 

% Check angular separation between sun vector S and slew plane 
alpha = pi/2 - acos(dot(S_G, e_G));         % coming out to 0 - check
% alpha = pi/4; 

%%%%%%
% IF angular separation is less than payload half-cone angle --> find phi2
% and phi3. Otherwise, slew is just phi1. 
%%%%%%

% unit_S = S_G/norm(S_G); 
S_PiPf_G = -cross(cross(S_G, e_G), e_G);    % sun projection vector G frame
S_PiPf_G = S_PiPf_G/norm(S_PiPf_G);         % sun projection --> unit vector 

% First slew around eigenaxis
phi1 = acos(dot(Pi_G, S_PiPf_G)) - ep;      % dot product of unit vectors 
if phi1 > pi/2 
    phi1_rem = phi1 - pi/2; 
end 

% Find P1 vector 
P1_P = [cos(phi1); sin(phi1); 0];           % P1 in P frame 
P1_G = G_DCM_P*P1_P;                        % P1 in G frame 

% Find P2 vector 
phiS = acos(dot(S_PiPf_G, Pi_G));           % angle btwn sun projection and Pi 
phiP2 = phiS + ep;                          % angle btwn P2 and Pi 
P2_P = [cos(phiP2); sin(phiP2); 0];         % P2 in P frame 
P2_G = G_DCM_P*P2_P;                        % P2 in G frame 

% Slew around sun vector via phi2 
if alpha == 0
    phi2 = pi; 
else 
    % Not sure how to get the one below: 
    theta = acos(dot(P1_G, S_G)); 
    phi2_M = 2*asin(sin(ep/2)/sin(theta/2)); 
    phi2_M2 = 2*asin(sin(ep)/(2*sin(theta/2)));     % REVISED FORMULA 

    % Junette's derived: 
    theta = acos(dot(S_G, P1_G));           % angle btwn sun and P1 vectors 
    P3_G = S_G*norm(P1_G)*cos(theta);       % P3 vector in G frame (S and P1 already unit vectors)
    P3P1_G = P1_G - P3_G;                   % vector from P3 to P1 
    P3P2_G = P2_G - P3_G;                   % vector from P3 to P2 
    phi2_P3 = acos(dot(P3P1_G/norm(P3P1_G), P3P2_G/norm(P3P2_G)));       % slew around sun vector 

    SP1_G = P1_G - S_G; 
    SP2_G = P2_G - S_G; 
    phi2_S = acos(dot(SP1_G/norm(SP1_G), SP2_G/norm(SP2_G))); 
    phi2 = phi2_S; 
end 

%%

% Slew around eigenvector via phi3 
% What if Pi and Pf overlap with P2 and P2? Need to ask Mohammad this ... 
phi3 = acos(dot(Pf_G, P2_G)); 
if phi3 > pi/2 
    phi3_rem = phi3 - pi/2; 
end 