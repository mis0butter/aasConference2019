% MAIN Example 
% Junette Hsin 
% Mohammad Ayoubi 

close all; 

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
G_DCM_N = eye(3);                           % DCM from the N to the G frame 
N_DCM_G = G_DCM_N';                         % G to N frame 
S_G = G_DCM_N*S_N;                          % Sun vector in G frame 
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

%% optional plotting routine to check things 

close all; 

plot_option = 0; 
if plot_option == 1
    figure()
        plot3([0 P1_G(1)], [0 P1_G(2)], [0 P1_G(3)], 'b'); 
        grid on; hold on
        plot(theta,sin(theta)./theta,'LineWidth',3) 
        plot3([0 P2_G(1)], [0 P2_G(2)], [0 P2_G(3)], 'b'); 
        plot3([0 S_G(1)], [0 S_G(2)], [0 S_G(3)], 'r'); 
        plot3([0 S_PiPf_G(1)], [0 S_PiPf_G(2)], [0 S_PiPf_G(3)], 'r'); 
        plot3([0 P3_G(1)], [0 P3_G(2)], [0 P3_G(3)], 'g'); 
        plot3([P3_G(1) P1_G(1)], [P3_G(2) P1_G(2)], [P3_G(3) P1_G(3)], 'g'); 
        plot3([P3_G(1) P2_G(1)], [P3_G(2) P2_G(2)], [P3_G(3) P2_G(3)], 'g'); 
        text(P3_G(1), P3_G(2), P3_G(3), ... 
            sprintf('    phi2 = %0.2f deg', phi2_P3*180/pi))
        title('Dot Product, Phi around P3') 
        xlabel('P perp_G') 
        ylabel('Pi_G')
        zlabel('e_G') 
        
    figure()
        plot3([0 P1_G(1)], [0 P1_G(2)], [0 P1_G(3)], 'b'); 
        grid on; hold on; 
        plot3([0 P2_G(1)], [0 P2_G(2)], [0 P2_G(3)], 'b'); 
        plot3([0 S_G(1)], [0 S_G(2)], [0 S_G(3)], 'g'); 
        plot3([0 S_PiPf_G(1)], [0 S_PiPf_G(2)], [0 S_PiPf_G(3)], 'r'); 
        plot3([S_G(1) P1_G(1)], [S_G(2) P1_G(2)], [S_G(3) P1_G(3)], 'g'); 
        plot3([S_G(1) P2_G(1)], [S_G(2) P2_G(2)], [S_G(3) P2_G(3)], 'g');
        plot3([P1_G(1) P2_G(1)], [P1_G(2) P2_G(2)], [P1_G(3) P2_G(3)], 'g');  
        text(S_G(1), S_G(2), S_G(3), ... 
            sprintf('    phi2_M = %0.2f deg \n     phi2_{M2} = %0.2f deg \n     phi2_S = %0.2f deg', ... 
            phi2_M*180/pi, phi2_M2*180/pi, phi2_S*180/pi))
        title('Chord Trigonometry, Phi around S') 
        xlabel('P perp_G') 
        ylabel('Pi_G')
        zlabel('e_G') 
end 

%%

% Slew around eigenvector via phi3 
% What if Pi and Pf overlap with P2 and P2? Need to ask Mohammad this ... 
phi3 = acos(dot(Pf_G, P2_G)); 
if phi3 > pi/2 
    phi3_rem = phi3 - pi/2; 
end 

%% Determine Phi 1 slew times

w0 = 0; 
t0 = 0; 
wf = 0; 

t1 = t0 + (wMax-w0)/aMax; 

% Mohammad's equation 
t2 = t1 - (1/wMax) * ... 
    ( phi1 - w0*(t1 - t0) - 0.5*aMax*(t1 - t0)^2 ... 
    - wMax*(wMax - wf)/aMax + (wMax - wf)^2/(2*aMax) );  

% % Phi1 times 
% t1_Phi1 = wMax/aMax; 
% t2_Phi1 = Phi1/wMax; 
% t3_Phi1 = t1_Phi1 + t2_Phi1; 

t3 = t2 - (wf - wMax)/aMax; 


%% Solve for attitue determination - first slew 
        
% t0 --> t1 
tEnd = t1 - t0; 
w_in = [    0;    0;      0]; 
q_in = [    0;      0;      0;      1];             % wrt G frame 
a = [ aMax;  0;  0];  
torque = inertia*a; 

[t1_phi1, y1_phi1] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [0, tEnd], [w_in; q_in]);

% t1 --> t2 
tEnd = t2 - t1; 
w_in = y1_phi1(end, 1:3)'; 
q_in = y1_phi1(end, 4:7)'; 
a = [0; 0; 0]; 
torque = inertia*a; 

[t2_phi1, y2_phi2] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [0, tEnd], [w_in; q_in]); 

% t2 --> t3 
tEnd = t3 - t2; 
w_in = y2_phi2(end, 1:3)'; 
q_in = y2_phi2(end, 4:7)'; 
a = [ -aMax; 0; 0]; 
torque = inertia*a; 

[t3_phi1, y3_phi1] = ode45(@(t, Z) gyrostat_cont(inertia, torque, Z), [0, tEnd], [w_in; q_in]); 

y_phi1 = [y1_phi1; y2_phi2(2:end, :); y3_phi1(2:end, :)]; 
w_phi1 = y_phi1(:, 1:3); 
q_phi1 = y_phi1(:, 4:7); 

ypr_phi1 = zeros(length(q_phi1), 3); 
for i = 1:max(size(q_phi1))
    ypr_phi1(i, :) = SpinCalc('QtoEA321', q_phi1(i, :), eps, 0); 
end 

t_phi1 = [t1_phi1; t1_phi1(end)+t2_phi1(2:end); t1_phi1(end)+t2_phi1(end)+t3_phi1(2:end)]; 

%% Plot phi1
% ylimits = get_ylimits(q); 
% ylimits = get_ylimits(w); 

% Plot 

if plot_option == 1
figure()
    plot(t_phi1, w_phi1) 
%     ylim(ylimits)
    legend('w1', 'w2', 'w3'); 
    ylabel('w (rad/s)') 
    xlabel('time (s)') 
    title('Angular Velocity Phi 1') 

figure()
    plot(t_phi1, q_phi1)
    legend('q1', 'q2', 'q3', 'q4'); 
%     ylim(ylimits)
    ylabel('quats') 
    xlabel('time (s)') 
    title('Quaternion Phi 1') 
    
figure()
    plot(t_phi1, ypr_phi1)
    legend('Yaw', 'Pitch', 'Roll'); 
    xlabel('time (s)') 
    ylabel('degrees') 
    title('Euler Angles Phi 1') 
end 
    
%% Determine Phi 2 slew times 

w0 = 0; 
t0 = 0; 
wf = 0; 

t1 = t0 + (wMax-w0)/aMax; 

% Mohammad's equation 
t2 = t1 - (1/wMax) * ... 
    ( phi2 - w0*(t1 - t0) - 0.5*aMax*(t1 - t0)^2 ... 
    - wMax*(wMax - wf)/aMax + (wMax - wf)^2/(2*aMax) );  

% % Phi1 times 
% t1_Phi1 = wMax/aMax; 
% t2_Phi1 = Phi1/wMax; 
% t3_Phi1 = t1_Phi1 + t2_Phi1; 

t3 = t2 - (wf - wMax)/aMax; 
    
%% Solve for attitude determination - second slew 

tEnd = t1 - t0; 
w_in = w_phi1(end, :)'; 
q_in = q_phi1(end, :)'; 
a = aMax*S_G; 
torque = inertia*a; 

% G --> N frame 
% N_Q_G = SpinCalc('DCMtoQ', N_DCM_G, eps, 1); 
% w_in = N_DCM_G*w_in; 
% q_in = N_Q_G*q_in; 
% torque_N = N_DCM_G*torque; 

%%
%%%%%%
% 
% What I need to happen: 
% for t0 --> t1: 
% 
% the direction of torque_G needs to be recalculated at every time step 

w = w_in; 
Q = q_in; 
dt = 1/16; 
int = 0; 
quat = zeros(t1/dt + 1, 4); 
ang_vel = zeros(t1/dt + 1, 3); 

for t = 0:dt:t1 
    nsteps = 10;
    for i = 1:nsteps
        dw = inv(inertia)*(torque - cross(w, inertia*w));
        dQ = 1/2*[Q(4), -Q(3), Q(2);
                  Q(3), Q(4), -Q(1);
                  -Q(2), Q(1), Q(4);
                  -Q(1), -Q(2), -Q(3)]*w;
        w = w + dw*dt/nsteps;
        Q = Q + dQ*dt/nsteps;
    end
    DCM = SpinCalc('QtoDCM', Q', eps, 0); 
    torque = DCM*torque; 
    inertia = DCM*inertia; 
    
    int = int + 1; 
    quat(int, :) = Q'; 
    ang_vel(int, :) = w'; 
    
end 
    
%%%%%%
%%

[t1_phi2, y1_phi2] = ode45(@(t,Z) gyrostat_cont(inertia, torque_N, Z), [0, tEnd], [w_in; q_in]);

tEnd = t2 - t1; 
w_in = y1_phi2(end, 1:3)'; 
q_in = y1_phi2(end, 4:7)'; 
a = [0; 0; 0]; 
torque = inertia*a; 

[t2_phi2, y2_phi2] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [0, tEnd], [w_in; q_in]);

tEnd = t3 - t1; 
w_in = y2_phi2(end, 1:3)'; 
q_in = y2_phi2(end ,4:7)'; 
a = -aMax*S_G; 
torque = inertia*a; 
torque_N = N_DCM_G*torque; 

[t3_phi2, y3_phi2] = ode45(@(t,Z) gyrostat_cont(inertia, torque_N, Z), [0, tEnd], [w_in; q_in]);

y_phi2 = [y1_phi2; y2_phi2(2:end, :); y3_phi2(2:end, :)]; 
w_phi2 = y_phi2(:, 1:3); 
q_phi2 = y_phi2(:, 4:7); 

ypr_phi2 = zeros(length(q_phi2), 3); 

for i = 1:max(size(q_phi2))
    ypr_phi2(i, :) = SpinCalc('QtoEA321', q_phi2(i, :), eps, 0); 
end 

t_phi2 = [t1_phi2; t1_phi2(end) + t2_phi2(2:end); t1_phi2(end) + t2_phi2(end) + t3_phi2(2:end)]; 

%% Plot phi2

w = y_phi2(:, 1:3); 
q = y_phi2(:, 4:7); 

if plot_option == 1
figure()
    plot(t_phi2, w)
%     ylim(ylimits)
    legend('w1', 'w2', 'w3'); 
    ylabel('w (rad/s)') 
    xlabel('time (s)') 
    title('Angular Velocity Phi 2') 

figure()
    plot(t_phi2, q)
    legend('q1', 'q2', 'q3', 'q4'); 
%     ylim(ylimits)
    ylabel('quats') 
    xlabel('time (s)') 
    title('Quaternion Phi 2') 
    
figure()
    plot(t_phi2, ypr_phi2)
    legend('Yaw', 'Pitch', 'Roll'); 
    xlabel('time (s)') 
    ylabel('degrees') 
    title('Euler Angles Phi 2') 
end 


%% Determine Phi 3 slew times

w0 = 0; 
t0 = 0; 
wf = 0; 

t1 = t0 + (wMax-w0)/aMax; 

% Mohammad's equation 
t2 = t1 - (1/wMax) * ... 
    ( phi1 - w0*(t1 - t0) - 0.5*aMax*(t1 - t0)^2 ... 
    - wMax*(wMax - wf)/aMax + (wMax - wf)^2/(2*aMax) );  

% % Phi1 times 
% t1_Phi1 = wMax/aMax; 
% t2_Phi1 = Phi1/wMax; 
% t3_Phi1 = t1_Phi1 + t2_Phi1; 

t3 = t2 - (wf - wMax)/aMax; 

%% Solve for attitue determination - third slew 
        
% t0 --> t1 
tEnd = t1 - t0; 
w_in = w_phi2(end, :)'; 
q_in = q_phi2(end, :)'; 
a = [ aMax;  0;  0];  
torque = inertia*a; 

[t1_phi3, y1_phi3] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [0, tEnd], [w_in; q_in]);

% t1 --> t2 
tEnd = t2 - t1; 
w_in = y1_phi3(end, 1:3)'; 
q_in = y1_phi3(end, 4:7)'; 
a = [0; 0; 0]; 
torque = inertia*a; 

[t2_phi3, y2_phi3] = ode45(@(t,Z) gyrostat_cont(inertia, torque, Z), [0, tEnd], [w_in; q_in]); 

% t2 --> t3 
tEnd = t3 - t2; 
w_in = y2_phi3(end, 1:3)'; 
q_in = y2_phi3(end, 4:7)'; 
a = [ -aMax; 0; 0]; 
torque = inertia*a; 

[t3_phi3, y3_phi3] = ode45(@(t, Z) gyrostat_cont(inertia, torque, Z), [0, tEnd], [w_in; q_in]); 

y_phi3 = [y1_phi3; y2_phi3(2:end, :); y3_phi3(2:end, :)]; 
w_phi3 = y_phi3(:, 1:3); 
q_phi3 = y_phi3(:, 4:7); 

ypr_phi3 = zeros(length(q_phi3), 3); 
for i = 1:max(size(q_phi3))
    ypr_phi3(i, :) = SpinCalc('QtoEA321', q_phi3(i, :), eps, 0); 
end 

t_phi3 = [t1_phi3; t1_phi3(end)+t2_phi3(2:end); t1_phi3(end)+t2_phi3(end)+t3_phi3(2:end)]; 

%% Plot phi3
% ylimits = get_ylimits(q); 
% ylimits = get_ylimits(w); 

% Plot 
if plot_option == 1
figure()
    plot(t_phi3, w_phi3) 
%     ylim(ylimits)
    legend('w1', 'w2', 'w3'); 
    ylabel('w (rad/s)') 
    xlabel('time (s)') 
    title('Angular Velocity Phi 3') 

figure()
    plot(t_phi3, q_phi3)
    legend('q1', 'q2', 'q3', 'q4'); 
%     ylim(ylimits)
    ylabel('quats') 
    xlabel('time (s)') 
    title('Quaternion Phi 3') 
    
figure()
    plot(t_phi3, ypr_phi3)
    legend('Yaw', 'Pitch', 'Roll'); 
    xlabel('time (s)') 
    ylabel('degrees') 
    title('Euler Angles Phi 3') 
end 

%% total slew stuff 

t_total = [ t_phi1; ... 
            t_phi1(end) + t_phi2(2:end); ... 
            t_phi1(end) + t_phi2(end) + t_phi3(2:end)]; 
        
w_total = [ w_phi1; ... 
            w_phi2(2:end, :); ... 
            w_phi3(2:end, :)]; 
        
q_total = [ q_phi1; ... 
            q_phi2(2:end, :); ... 
            q_phi3(2:end, :)]; 
        
ypr_total = zeros(length(q_total), 3); 
for i = 1:max(size(q_total))
    ypr_total(i, :) = SpinCalc('QtoEA321', q_total(i, :), eps, 0); 
end 
        
%% plot total stuff 
        
plot_total = 1; 
if plot_total == 1

figure()
plot(t_total, w_total) 
%     ylim(ylimits)
    grid on
    legend('w1', 'w2', 'w3'); 
    ylabel('w (rad/s)') 
    xlabel('time (s)') 
    title('Angular Velocity Total') 

figure()
    plot(t_total, q_total)
    grid on 
    legend('q1', 'q2', 'q3', 'q4'); 
%     ylim(ylimits)
    ylabel('quats') 
    xlabel('time (s)') 
    title('Quaternion Total') 
    
figure()
    plot(t_total, ypr_total)
    grid on 
    legend('Yaw', 'Pitch', 'Roll'); 
    xlabel('time (s)') 
    ylabel('degrees') 
    title('Euler Angles Total') 
    
end 

%% 
% 
% function ylimits = get_ylimits(data) 
% % Set y limits of axis 
% 
% rng = range(data(:)); 
% midp = min(data(:)) + rng/2; 
% ylimits = [ midp - 1.2*rng/2, midp + 1.2*rng/2 ]; 
% 
% end 