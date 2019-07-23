function [t1, t2, t3] = find_slew_times(a_max, w_max, phi, w0, wf, t0)

% find slew times 
% Mohammad Ayoubi 
% Junette Hsin 
    
%% Determine Phi slew times - presentation / paper equations 

t1 = t0 + (w_max-w0)/a_max; 

t2 = t1 + (1/w_max) * ... 
    ( phi - w0*(t1 - t0) - 0.5*a_max*(t1 - t0)^2 ... 
    - w_max*(w_max - wf)/a_max + (w_max - wf)^2/(2*a_max) );  

t3 = t2 - (wf - w_max)/a_max; 

%% Determine Phi slew times - Mohammad's scratch derived equations 
% 
% t1 = t0 + (w_max-w0)/a_max; 
% 
% t2 = t1 + (1/w_max) * ... 
%     ( phi + w_max*(wf - w_max)/a_max + 0.5*a_max*((wf - w_max)/a_max)^2 ... 
%     - a_max/2*(t1 - t0)^2 - w0*(t1 - t0) ); 
% 
% t3 = t2 - (wf - w_max)/a_max; 

end 