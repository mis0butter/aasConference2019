% torque envelope 

% alpha = 60 deg 
a = 60*pi/180; 

% Wheel axes from SC axes transformation matrix 
whl_mx_sc = [ ... 
    sqrt(2)/2 * sin(a)      -sqrt(2)/2 * sin(a)     cos(a)  ; 
    sqrt(2)/2 * sin(a)      sqrt(2)/2 * sin(a)      cos(a)  ; 
    -sqrt(2)/2 * sin(a)     -sqrt(2)/2 * sin(a)     cos(a)  ; 
    -sqrt(2)/2 * sin(a)     sqrt(2)/2 * sin(a)      cos(a)  ]; 

% SC axes from wheel transformation matrix 
sc_mx_whl = (whl_mx_sc' * whl_mx_sc) ^-1 * whl_mx_sc'; 

% SC inertia 
inertia_SC = [ 100     0       0; 
               0       100     0; 
               0       0       100]; 
inertia_whl = [ 4e-4    0       0;
                0       4e-4    0; 
                0       0       4e-4]; 
           
% four wheel torque capability --> torque ON spacecraft 
t = -1 : 0.1 : 1; 



% iterate through each wheel, calculate torque 
count = 0; 
for i = 1:numel(t)
    for j = 1:numel(t) 
        for k = 1:numel(t) 
            for m = 1:numel(t) 
                
                count = count + 1; 
                
                % Available torque 
                T_1234(:,count) = t(i) * sc_mx_whl(:,1) + t(j) * sc_mx_whl(:,2) ... 
                    + t(k) * sc_mx_whl(:,3) + t(m) * sc_mx_whl(:,4); 
                
                Tnorm_1234(count) = norm(T_1234(:,count)); 
                ijk_1234 = [count, i, j, k]; 
                
                % Available acceleration 
                wdot(:,count) = inv(inertia_SC) * T_1234(:, count); 
                wdot_norm(count) = norm(wdot(:,count)); 
                
            end
        end
    end
end

