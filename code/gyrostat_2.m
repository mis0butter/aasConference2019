function f=gyrostat_2(kp,kd,x,t,omega_c,omegadot_c,q_c,Jnom)
global torq
% Kinematic and EoMs of a Rigid Gyrostat

q=[x(1);x(2);x(3);x(4)]; % Gyrosdtat quaternions: [q1,q2,q3,q4]
omega=[x(5);x(6);x(7)]; % Gyrostat angular velocity components
wh=[x(8);x(9);x(10)]; % Wheel momenta
q_cc=[x(11);x(12);x(13);x(14)]; % computed commanded quarternions

torq=NLC_fun(kp,kd,omega,omega_c,omegadot_c,q,q_cc,Jnom);
% Integration Functions
skew_omega=[  0      -omega(3)  omega(2);...
            omega(3)    0       -omega(1);...
            -omega(2)  omega(1)     0];
Omega=[-skew_omega omega;-omega' 0]; % Eq. (2.97a)
skew_omega_c=[  0      -omega_c(3)  omega_c(2);...
            omega_c(3)    0       -omega_c(1);...
            -omega_c(2)  omega_c(1)     0];
Omega_c=[-skew_omega_c omega_c;-omega_c' 0]; % Eq. (2.97a)
f(1:4,:)=0.5*Omega*q;  % Eq. (7.1)
f(5:7,:)=Jnom\(-skew_omega*Jnom*omega+torq);% Eq. (7.16a)
f(8:10,:)=-skew_omega*wh-torq;% Wheel torque; Eq. (7.16b)
% Computing the commanded quarternions from the commanded angular velocity
f(11:14,:)=0.5*Omega_c*q_cc;% commanded quarternion

