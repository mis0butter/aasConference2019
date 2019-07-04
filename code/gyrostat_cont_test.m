% gyrostat_cont_test function 

w_in = [0.5; 0; 0]; 
q_in = [0; 0; 0; 1]; 
torque_x = [1; 0; 0]; 
inertia = [ 408     0       0; 
            0       427     0; 
            0       0       305]; 
        
[t, y] = ode45(@(t,Z) gyrostat_cont(inertia, torque_x, Z), [0, 100], [w_in; q_in]);
