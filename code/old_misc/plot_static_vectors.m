% plot static initial, final, sun vectors 
    
plot3([0 Pi_G0(1)], [0 Pi_G0(2)], [0 Pi_G0(3)], 'b'); 
plot3([0 Pf_G0(1)], [0 Pf_G0(2)], [0 Pf_G0(3)], 'b'); 
plot3([0 P1_G0(1)], [0 P1_G0(2)], [0 P1_G0(3)], 'b'); 
plot3([0 P2_G0(1)], [0 P2_G0(2)], [0 P2_G0(3)], 'b'); 
plot3([0 e_G0(1)], [0 e_G0(2)], [0 e_G0(3)], 'b'); 

plot(theta,sin(theta)./theta,'LineWidth',3) 

plot3([0 S_G0(1)], [0 S_G0(2)], [0 S_G0(3)], 'r'); 
plot3([0 S_PiPf_G0(1)], [0 S_PiPf_G0(2)], [0 S_PiPf_G0(3)], 'r'); 
plot3([0 P3_G0(1)], [0 P3_G0(2)], [0 P3_G0(3)], 'g'); 

plot3([P3_G0(1) P1_G0(1)], [P3_G0(2) P1_G0(2)], [P3_G0(3) P1_G0(3)], 'g'); 
plot3([P3_G0(1) P2_G0(1)], [P3_G0(2) P2_G0(2)], [P3_G0(3) P2_G0(3)], 'g'); 