function B_q_C = quat_diff(A_q_B, A_q_C)

qc = A_q_C; 
qb = A_q_B; 

qb_skew = [  qb(4)  qb(3) -qb(2) -qb(1); ...
            -qb(3)  qb(4)  qb(1) -qb(2); ...
             qb(2) -qb(1)  qb(4) -qb(3); ... 
             qb(1)  qb(2)  qb(3)  qb(4)]; 
B_q_C = qb_skew * qc; 
