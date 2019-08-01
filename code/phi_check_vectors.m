% Script to check phi 1 vectors 

figure()
    
    plot3([0 Pi_G0(1)], [0 Pi_G0(2)], [0 Pi_G0(3)], 'b'); 
    hold on; grid on 
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
    
    text(P3_G0(1), P3_G0(2), P3_G0(3), ... 
        sprintf('    phi2 = %0.2f deg', phi2_P3*180/pi))
    text(Pi_G0(1), Pi_G0(2), Pi_G0(3), sprintf(' Pi')) 
    text(Pf_G0(1), Pf_G0(2), Pf_G0(3), sprintf(' Pf')) 
    text(e_G0(1), e_G0(2), e_G0(3), sprintf(' e')) 
    text(P1_G0(1), P1_G0(2), P1_G0(3), sprintf(' P1')) 
    text(P2_G0(1), P2_G0(2), P2_G0(3), sprintf(' P2')) 
    text(S_G0(1), S_G0(2), S_G0(3), sprintf(' sun')) 
    text(S_PiPf_G0(1), S_PiPf_G0(2), S_PiPf_G0(3), sprintf(' sun proj')) 
    title('Dot Product, Phi around P3') 
    xlabel('G0_x')
    ylabel('G0_y') 
    zlabel('G0_z') 
%     xlabel('P perp_G0') 
%     ylabel('Pi_G0')
%     zlabel('e_G0') 

figure()
    plot3([0 P1_G0(1)], [0 P1_G0(2)], [0 P1_G0(3)], 'b'); 
    grid on; hold on; 
    plot3([0 P2_G0(1)], [0 P2_G0(2)], [0 P2_G0(3)], 'b'); 
    plot3([0 Pi_G0(1)], [0 Pi_G0(2)], [0 Pi_G0(3)], 'b'); 
    plot3([0 Pf_G0(1)], [0 Pf_G0(2)], [0 Pf_G0(3)], 'b'); 
    plot3([0 e_G0(1)], [0 e_G0(2)], [0 e_G0(3)], 'b'); 
    plot3([0 S_G0(1)], [0 S_G0(2)], [0 S_G0(3)], 'g'); 
    plot3([0 S_PiPf_G0(1)], [0 S_PiPf_G0(2)], [0 S_PiPf_G0(3)], 'r'); 
    plot3([S_G0(1) P1_G0(1)], [S_G0(2) P1_G0(2)], [S_G0(3) P1_G0(3)], 'g'); 
    plot3([S_G0(1) P2_G0(1)], [S_G0(2) P2_G0(2)], [S_G0(3) P2_G0(3)], 'g');
    plot3([P1_G0(1) P2_G0(1)], [P1_G0(2) P2_G0(2)], [P1_G0(3) P2_G0(3)], 'g');  
    
    text(S_G0(1), S_G0(2), S_G0(3), ... 
        sprintf('    phi2_M = %0.2f deg \n     phi2_{M2} = %0.2f deg \n     phi2_S = %0.2f deg', ... 
        phi2_M*180/pi, phi2_M2*180/pi, phi2_S*180/pi))
    text(Pi_G0(1), Pi_G0(2), Pi_G0(3), sprintf(' Pi')) 
    text(Pf_G0(1), Pf_G0(2), Pf_G0(3), sprintf(' Pf')) 
    text(e_G0(1), e_G0(2), e_G0(3), sprintf(' e'))  
    text(P1_G0(1), P1_G0(2), P1_G0(3), sprintf(' P1')) 
    text(P2_G0(1), P2_G0(2), P2_G0(3), sprintf(' P2')) 
    text(S_G0(1), S_G0(2), S_G0(3), sprintf(' sun'))  
    text(S_PiPf_G0(1), S_PiPf_G0(2), S_PiPf_G0(3), sprintf(' sun proj')) 
    title('Chord Trigonometry, Phi around S') 
    xlabel('G0_x')
    ylabel('G0_y') 
    zlabel('G0_z') 
%     xlabel('P perp_G0') 
%     ylabel('Pi_G0')
%     zlabel('e_G0') 
        



%     plot3([0 P1_G0(1)], [0 P1_G0(2)], [0 P1_G0(3)], 'b'); 
%     grid on; hold on
%     plot(theta,sin(theta)./theta,'LineWidth',3) 
%     plot3([0 P2_G0(1)], [0 P2_G0(2)], [0 P2_G0(3)], 'b'); 
%     plot3([0 Pi_G0(1)], [0 Pi_G0(2)], [0 Pi_G0(3)], 'b'); 
%     plot3([0 Pf_G0(1)], [0 Pf_G0(2)], [0 Pf_G0(3)], 'b'); 
%     plot3([0 e_G0(1)], [0 e_G0(2)], [0 e_G0(3)], 'b'); 
%     plot3([0 S_G0(1)], [0 S_G0(2)], [0 S_G0(3)], 'r'); 
%     plot3([0 S_PiPf_G0(1)], [0 S_PiPf_G0(2)], [0 S_PiPf_G0(3)], 'r'); 
%     plot3([0 P3_G0(1)], [0 P3_G0(2)], [0 P3_G0(3)], 'g'); 
%     plot3([P3_G0(1) P1_G0(1)], [P3_G0(2) P1_G0(2)], [P3_G0(3) P1_G0(3)], 'g'); 
%     plot3([P3_G0(1) P2_G0(1)], [P3_G0(2) P2_G0(2)], [P3_G0(3) P2_G0(3)], 'g'); 