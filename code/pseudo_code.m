% MAIN Example 
close all; 

% Non-changing parameters; should be function inputs 
inertia = [ 408     0       0; 
            0       427     0; 
            0       0       305]; 
Pi = [1; 0; 0]; 					% Pi = unit vector of initial point in the G frame 
Pf = [0; 1; 0];  					% Pf = unit vector of the final point in the G frame 
S_N = [cosd(45); sind(45); 0];      % S = unit vector of sun vector in the N frame 
G_DCM_N = eye(3);                   % Transformation from the N to the G frame 
N_DCM_G = G_DCM_N';                 % G to N frame 
S = G_DCM_N*S_N;                    % Sun vector in G frame 
ep = pi/12;                         % payload half-cone angle. pi/12 rad = 15 deg  

% Calculate normal vector of slew plane 
e = cross(Pi, Pf) / norm(cross(Pi, Pf));    % eigenaxis of PiPf plane, in G frame 

% Check angular separation between sun vector S and slew plane 
alpha = pi/2 - acos(dot(S, e)); 

% If angular separation is less than payload half-cone angle
if alpha < ep 	
    unit_S = S/norm(S); 
    S_PiPf = cross(cross(unit_S, e), e);    % projection of sun onto slew plane, G frame 
end 

% First slew around eigenaxis
phi_1 = acos(dot(Pi, S_PiPf)) - ep;       % dot product of unit vectors 

% Slew around sun vector via phi_2 
if alpha == 0
    phi_2 = pi; 
else 
    theta = acos(dot(Pi, S_N)); 
    phi_2 = 2*asin(sin(ep)/sin(theta)); 
end 

% Slew around eigenaxis through phi_3 
% if dot(Pf, 
%     
% else 
%     
% end 

%%
        
t = 100; 

phi_1_accel = [ 2*(phi_1/2)/t^2;  0;  0];  
torque = inertia*phi_1_accel; 

w_in = [    0.5;    0;      0]; 
q_in = [    0;      0;      0;      1]; 

[t, y1] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [0, 100], [w_in; q_in]);

w_in = y1(end, 1:3)'; 
q_in = y1(end, 4:7)'; 

[t, y2] = ode45(@(t,Z) gyrostat_cont(inertia, -torque, Z), [0, 100], [w_in; q_in]); 

y = [y1; y2]; 

% Plot 
figure()
plot(y(:, 1:3)) 
legend('w1', 'w2', 'w3'); 

figure()
plot(y(:, 4:7))
legend('q1', 'q2', 'q3', 'q4'); 