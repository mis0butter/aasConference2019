function DCM = quat2DCM(q)
% Junette Hsin 
% From Keith Reckdahl's memo 

% Describing the orientation of B in A --> A_DCM_B

C11 = 1 - 2*q(2)^2 - 2*q(3)^2; 
C12 = 2*(q(1)*q(2) - q(3)*q(4)); 
C13 = 2*(q(1)*q(3) + q(2)*q(4)); 

C21 = 2*(q(1)*q(2) + q(3)*q(4)); 
C22 = 1 - 2*q(1)^2 - 2*q(3)^2; 
C23 = 2*(q(2)*q(3) - q(1)*q(4)); 

C31 = 2*(q(1)*q(3) - q(2)*q(4)); 
C32 = 2*(q(1)*q(4) + q(2)*q(3)); 
C33 = 1 - 2*q(1)^2 - 2*q(2)^2; 

DCM = [C11 C12 C13; C21 C22 C23; C31 C32 C33]; 


