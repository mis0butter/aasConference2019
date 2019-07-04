

% Simulation parameters 
sim_rate = 0.01;			% sim rate = 100 Hz 
sim_duration = 60*60;  			% sim duration = 60 minutes = 3600 seconds 
t =  0:sim_rate:sim_duration; 		% defining time vector  

% Non-changing parameters 
Pi = [1; 0; 0]; 					% Pi = unit vector of initial point in the G frame 
Pf = [0; 1; 0];  					% Pf = unit vector of the final point in the G frame 
S = [1; 0; 0];                      % S = unit vector of sun vector in the N frame 

G_N_DCM = eye(3); 

ep = pi/12; 					% payload half-cone angle. pi/12 rad = 15 deg  
e = cross(Pi, Pf) / norm(cross(Pi, Pf)); 			% eigenaxis of PiPf plane 

% Check angular separation between sun vector S and the plane of P 
alpha = pi/2 - acos(dot(S, e)); 

% If angular separation is less than payload half-cone angle
if alpha < ep 	
        S_PiPf = S*cos(alpha); 
end 

phi_1a = acos(dot(Pi, S_PiPf)) - ep; 
phi_1b = acos(dot(Pi, S_PiPf)) - ep - 2*pi; 

% each iteration of i represents simulation running at sim rate 
for i = 1:max(size(t))/3
	
     phi_1 = acosd(dot(Pi, S_PiPf)) - ep; 
           
		% Slew around the Sun vector 
		if alpha == 0 
			phi_2 = pi; 
		else 
			phi_2 = 2*atan(ep/sin(alpha)); 
		end 

	% Else, if angular separation greater than payload half-cone angle 
	else 
% 		% Slew around the eigenaxis. If before sun slew: 
% 		if before sun slew: 
% 			phi_1 =  
% 		if after sun slew: 
% 			phi_3 =  
	end 

	% Perform slew around phi at sim_rate 
% 	slew ... 

end 
