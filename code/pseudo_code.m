% MAIN Example 
close all; 

% Non-changing parameters; should be function inputs 
inertia = [ 408     0       0; 
            0       427     0; 
            0       0       305]; 
Pi_G = [1; 0; 0];                           % Pi = unit vector of initial point in the G frame 
Pf_G = [0; 1; 0];                           % Pf = unit vector of the final point in the G frame 
S_N = [cosd(45); sind(45); 0];              % S = unit vector of sun vector in the N frame 
G_DCM_N = eye(3);                           % DCM from the N to the G frame 
N_DCM_G = G_DCM_N';                         % G to N frame 
S_G = G_DCM_N*S_N;                          % Sun vector in G frame 
ep = pi/12;                                 % payload half-cone angle. pi/12 rad = 15 deg  

% Calculate normal vector of slew plane 
e_G = cross(Pi_G, Pf_G) / norm(cross(Pi_G, Pf_G));  % eigenaxis of PiPf plane, in G frame 
Pperp_G = cross(e_G, Pi_G/norm(Pi_G));              % perpendicular vector to e and Pi, slew plane, G frame 
P_DCM_G = [(Pi_G/norm(Pi_G))'; Pperp_G'; e_G'];     % DCM from G frame to P (slew plane) frame 
G_DCM_P = P_DCM_G';                         % P to G frame 

% Check angular separation between sun vector S and slew plane 
alpha = pi/2 - acos(dot(S_G, e_G)); 

%%%%%%
% IF angular separation is less than payload half-cone angle --> find phi2
% and phi3. Otherwise, slew is just phi1. 
%%%%%%

unit_S = S_G/norm(S_G); 
S_PiPf_G = cross(cross(unit_S, e_G), e_G);  % sun projection vector G frame 

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
%     theta = acos(dot(Pi_G, S_N)); 
%     phi2 = 2*asin(sin(ep)/sin(theta)); 

    % Junette's derived: 
    theta = acos(dot(S_G, P1_G));           % angle btwn sun and P1 vectors 
    P3_G = S_G*P1_G*cos(theta);             % P3 vector in G frame (S and P1 already unit vectors)
    P3P1_G = P1_G - P3_G;                   % vector from P3 to P1 
    P3P2_G = P2_G - P3_G;                   % vector from P3 to P2 
    phi2 = acos(dot(P3P1_G, P3P2_G));       % slew around sun vector 
end 

% Slew around eigenvector via phi3 
% What if Pi and Pf overlap with P2 and P2? Need to ask Mohammad this ... 
phi3 = acos(dot(Pf_G, P2_G)); 
if phi3 > pi/2 
    phi3_rem = phi3 - pi/2; 
end 



%%
        
t = 100; 

phi_1_accel = [ 2*(phi1/2)/t^2;  0;  0];  
torque = inertia*phi_1_accel; 

w_in = [    0.5;    0;      0]; 
q_in = [    0;      0;      0;      1]; 

[t1, y1] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [0, 100], [w_in; q_in]);

w_in = y1(end, 1:3)'; 
q_in = y1(end, 4:7)'; 

[t2, y2] = ode45(@(t,Z) gyrostat_cont(inertia, -torque, Z), [0, 100], [w_in; q_in]); 

y = [y1; y2]; 
t = [t1; t2]; 

%%

w = y(:, 1:3); 
ylimits = get_ylimits(w); 

% Plot 
figure()
    plot(t, w) 
    ylim(ylimits)
    legend('w1', 'w2', 'w3'); 
    ylabel('w (rad/s)') 
    xlabel('time (s)') 
    title('Angular Velocity') 

q = y(:, 4:7); 
ylimits = get_ylimits(q); 

figure()
    plot(t, q)
    legend('q1', 'q2', 'q3', 'q4'); 
    ylim(ylimits)
    ylabel('quats') 
    xlabel('time (s)') 
    title('Quaternions') 

%% 
function ylimits = get_ylimits(data) 
% Set y limits of axis 

rng = range(data(:)); 
midp = min(data(:)) + rng/2; 
ylimits = [ midp - 1.2*rng/2, midp + 1.2*rng/2 ]; 

end 