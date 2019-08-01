

%% Inputs 

% Non-changing parameters; should be function inputs 
inertia_SC = [ 408     0       0; 
               0       427     0; 
               0       0       305]; 
inertia_w = [  20   0    0; 
               0    20   0; 
               0    0    20  ]; 
aMax = 1;                                  % Maximum acceleration, rad/s^2
wMax = 1;                                  % Maximum angular velocity, rad/s
ep = pi/12;                                 % payload half-cone angle. pi/12 rad = 15 deg  

% Initial points / vectors 
Pi_G0 = [rand(1); rand(1); rand(1)];                           % Pi = unit vector of initial point in the G frame 
Pi_G0 = Pi_G0/norm(Pi_G0); 
Pf_G0 = [rand(1); rand(1); rand(1)];                           % Pf = unit vector of the final point in the G frame 
Pf_G0 = Pf_G0/norm(Pf_G0); 
while acos(dot(Pi_G0, Pf_G0)) < ep*2
    Pf_G0 = [rand(1); rand(1); rand(1)]; 
    Pf_G0 = Pf_G0/norm(Pf_G0); 
end 

% Calculate normal vector of slew plane 
e_G0 = cross(Pf_G0, Pi_G0) / norm(cross(Pf_G0, Pi_G0));  % eigenaxis of PiPf plane, in G frame 
Pperp_G0 = cross(Pi_G0, e_G0);              % perpendicular vector to e and Pi, slew plane, G frame 
P_DCM_G0 = [Pi_G0'; Pperp_G0'; e_G0'];     % DCM from G frame to P (slew plane) frame 
G0_DCM_P = P_DCM_G0';                         % P to G frame 
% Defining inertial frames 
% G0_DCM_N = angle2dcm(1, 1, 1);         % DCM from the N to the G frame 
G0_DCM_N = eye(3); 
N_DCM_G0 = G0_DCM_N';                         % G to N frame - initial!!! G frame will change throughout sim 

% Sun vector stuff 
S_N = [cosd(75); sind(75); 1];              % S = unit vector of sun vector in the N frame 
S_N = S_N/norm(S_N);                        % normalizing sun vector 
S_G0 = G0_DCM_N*S_N;                          % Sun vector in G frame 

% Check angular separation between sun vector S and slew plane 
alpha = pi/2 - acos(dot(S_G0, e_G0));         % coming out to 0 - check

%%%%%%
% IF angular separation is less than payload half-cone angle --> while loop
% until alpha < ep. for simulation!!! 
%%%%%%

while alpha*180/pi > ep 
    disp('alpha > ep') 
    S_N = [rand(1); rand(1); rand(1)];  
    S_G0 = G0_DCM_N*S_N;               
    alpha = pi/2 - acos(dot(S_G0, e_G0));         % coming out to 0 - check           
end 

% FINALLY - Sun projection onto G frame 
S_PiPf_G0 = -cross(cross(S_G0, e_G0), e_G0);    % sun projection vector G frame
S_PiPf_G0 = S_PiPf_G0/norm(S_PiPf_G0);         % sun projection --> unit vector 

%% Calculate slew angles 

% Calculate the treshold triangle 
% phi_tt = 2*wMax*sqrt(aMax^2 + wMax^2);      % Threshold triangle!!! 
phi_tt = wMax^2/aMax; 

% First slew around eigenaxis
phi1 = acos(dot(Pi_G0, S_PiPf_G0)) - ep;      % dot product of unit vectors 
if phi1 > pi/2 
    phi1_rem = phi1 - pi/2; 
end 

% Find P1 vector in G0 frame 
P1_P = [cos(phi1); sin(phi1); 0];           % P1 in P frame 
P1_G0 = G0_DCM_P*P1_P;                        % P1 in G frame 

% Find P2 vector 
phiS = acos(dot(S_PiPf_G0, Pi_G0));           % angle btwn sun projection and Pi 
phiP2 = phiS + ep;                          % angle btwn P2 and Pi 
P2_P = [cos(phiP2); sin(phiP2); 0];         % P2 in P frame 
P2_G0 = G0_DCM_P*P2_P;                        % P2 in G frame 

% Slew around sun vector via phi2 
if alpha == 0
    phi2 = pi; 
else 
    % Not sure how to get the one below: 
    theta = acos(dot(P1_G0, S_G0)); 
    phi2_M = 2*asin(sin(ep/2)/sin(theta/2)); 
    phi2_M2 = 2*asin(sin(ep)/(2*sin(theta/2)));     % REVISED FORMULA 

    % Junette's derived: 
    theta = acos(dot(S_G0, P1_G0));           % angle btwn sun and P1 vectors 
    P3_G0 = S_G0*norm(P1_G0)*cos(theta);       % P3 vector in G frame (S and P1 already unit vectors)
    P3P1_G0 = P1_G0 - P3_G0;                   % vector from P3 to P1 
    P3P2_G0 = P2_G0 - P3_G0;                   % vector from P3 to P2 
    phi2_P3 = acos(dot(P3P1_G0/norm(P3P1_G0), P3P2_G0/norm(P3P2_G0)));       % slew around sun vector 

    SP1_G0 = P1_G0 - S_G0; 
    SP2_G0 = P2_G0 - S_G0; 
    phi2_S = acos(dot(SP1_G0/norm(SP1_G0), SP2_G0/norm(SP2_G0))); 
    phi2 = phi2_S; 
end 

%%

% Slew around eigenvector via phi3 
% What if Pi and Pf overlap with P2 and P2? Need to ask Mohammad this ... 
phi3 = acos(dot(Pf_G0, P2_G0)); 
if phi3 > pi/2 
    phi3_rem = phi3 - pi/2; 
end 