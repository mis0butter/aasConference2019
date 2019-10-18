% Junette Hsin 
% 2019 August 9 
% Inputs 

%% Non-changing parameters; should be function inputs 

inertia_SC = [ 100     0       0; 
               0       100     0; 
               0       0       100]; 
inertia_w = [  20   0    0; 
               0    20   0; 
               0    0    20  ]; 

ep = pi/12;                                 % payload half-cone angle. pi/12 rad = 15 deg  
% aMax = 0.2;                                  % Maximum acceleration, rad/s^2
% wMax = 0.1;                                  % Maximum angular velocity, rad/s
aMax = 10*rand; 
wMax = 10*rand; 

% Ensure discretization of slew times is greater than 10 ms
while wMax/aMax < 0.1
    aMax = 10*rand; 
    wMax = 10*rand; 
end 

%% Create random points 

% Initial points / vectors 
Pi_N = [rand*(-1)^round(rand); rand*(-1)^round(rand); rand*(-1)^round(rand)]; Pi_N = Pi_N / norm(Pi_N);          
Pf_N = [rand*(-1)^round(rand); rand*(-1)^round(rand); rand*(-1)^round(rand)]; Pf_N = Pf_N / norm(Pf_N);  
while acos(dot(Pi_N, Pf_N)) < pi/3
    Pf_N = [rand*(-1)^round(rand); rand*(-1)^round(rand); rand*(-1)^round(rand)]; 
    Pf_N = Pf_N / norm(Pf_N); 
end 

% Calculate normal vector of slew plane 
e_N = cross(Pf_N, Pi_N) / norm(cross(Pf_N, Pi_N));  % eigenaxis of PiPf plane, in N frame 
Pperp_N = cross(Pi_N, e_N);                         % perpendicular vector to e and Pi, slew plane, N frame 
P_DCM_N = [Pi_N'; Pperp_N'; e_N'];                  % DCM from N frame to P (slew plane) frame 
N_DCM_P = P_DCM_N';                                 % P to N frame 

% Create sun vector 
[alpha, theta_Pi_Sproj, theta_Sproj_Pf, theta_Pi_Pf, S_N, S_PiPf_N] = ... 
    sun_vector(e_N, Pi_N, Pf_N); 

% IF angular separation is less than payload half-cone angle --> while loop
% until alpha < ep. for simulation!!! 
while abs(alpha) > ep  || theta_Sproj_Pf < ep || theta_Pi_Sproj < ep || ... 
        theta_Pi_Sproj > theta_Pi_Pf || theta_Sproj_Pf > theta_Pi_Pf || theta_Sproj_Pf + theta_Pi_Sproj > pi

    [alpha, theta_Pi_Sproj, theta_Sproj_Pf, theta_Pi_Pf, S_N, S_PiPf_N] = ... 
        sun_vector(e_N, Pi_N, Pf_N); 
    
end 

%% Calculate slew angles 

% Calculate the threshold angle 
phi_t = wMax^2/aMax; 

% Find phi1 
phi1 = acos(dot(Pi_N, S_PiPf_N)) - ep;      % dot product of unit vectors 

% Find P1 vector in N frame 
P1_P = [ cos(phi1); sin(phi1); 0 ];         % P1 in P frame 
P1_N = N_DCM_P*P1_P;                        % P1 in N frame 

% Find P2 vector 
phiS = acos(dot(S_PiPf_N, Pi_N));           % angle btwn sun projection and Pi 
phiP2 = phiS + ep;                          % angle btwn P2 and Pi 
P2_P = [cos(phiP2); sin(phiP2); 0];         % P2 in P frame 
P2_N = N_DCM_P*P2_P;                        % P2 in N frame 

%% FIND PHI2

% Slew around sun vector via phi2 

if alpha == 0
    
    phi2 = pi; 
    theta = acos(dot(S_N, P1_N));          % angle btwn sun and P1 vectors 
    P3_N = S_N * norm(P1_N)*cos(theta);     % P3 vector in N frame (S and P1 already unit vectors)
    
else 
    
    % Not sure how to get the one below: 
    theta = acos(dot(P1_N, S_N)); 
    phi2_M = 2*asin(sin(ep/2)/sin(theta/2)); 
    phi2_M2 = 2*asin(sin(ep)/(2*sin(theta/2)));     % REVISED FORMULA 

    % Junette's derived: 
    theta = acos(dot(S_N, P1_N));           % angle btwn sun and P1 vectors 
    P3_N = S_N*norm(P1_N)*cos(theta);       % P3 vector in G frame (S and P1 already unit vectors)
    P3P1_N = P1_N - P3_N;                   % vector from P3 to P1 
    P3P2_N = P2_N - P3_N;                   % vector from P3 to P2 
    phi2_P3 = acos(dot(P3P1_N/norm(P3P1_N), P3P2_N/norm(P3P2_N)));       % slew around sun vector 

    SP1_N = P1_N - S_N; 
    SP2_N = P2_N - S_N; 
    phi2_S = acos(dot(SP1_N/norm(SP1_N), SP2_N/norm(SP2_N))); 
    
    % Another phi2
    top = dot(S_N, cross(P1_N, P2_N)); 
    bot = dot(P1_N, P2_N) - dot(S_N, P1_N)*dot(S_N, P2_N);
    phi2_M3 = abs( atan2 ( top, bot ) ); 
    
    % Last one? This doesn't work
    theta = acos(dot(P1_N, S_N)); 
    top = (pi/2 - alpha)*sin(ep); 
    bot = cos(ep) - cos(theta)*cos(alpha); 
    phi2_M4 = 2*abs(atan2(top, bot)); 

    % The chosen phi2! 
    phi2 = phi2_M2; 

end 


%% Find phi3

% Slew around eigenvector via phi3 
phi3 = acos(dot(Pf_N, P2_N)); 
if phi3 > pi/2 
    phi3_rem = phi3 - pi/2; 
end 
