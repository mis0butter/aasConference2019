t_0=0;

phidot_t0=0 ; %need to update later.  start with a slew from rest

phi_tf=phi1*100;

phidot_tf=0;

phidot_max = 0.1;

phiddot_max = 0.001;  %need to get this number from DG work

%switching times

T_1 = t_0 + (phidot_max-phidot_t0)/phiddot_max;

T_2 = T_1 + 1/phidot_max * ( phi_tf - phidot_t0*(T_1-t_0) - 1/2*phiddot_max*(T_1-t_0)^2 - phidot_max*(phidot_max-phidot_tf)/phiddot_max + ((phidot_max - phidot_tf)^2)/2/phiddot_max );

T_F = T_1 + 1/phidot_max * (phi_tf - phidot_t0*(T_1-t_0)- 1/2*phiddot_max*(T_1-t_0)^2 + ((phidot_max - phidot_tf)^2)/2/phiddot_max);

 

 

 

time(1)=0;

ndt_1 = ceil(T_F/0.0625);

for i=2:ndt_1+1,

    time(i) = time(i-1)+0.0625;

    if (time(i) >= t_0 && time(i) <= T_1)

        phiddot_t(i)=phiddot_max;

        phidot_t(i) = phidot_t0 + phiddot_max*(time(i) - t_0);

        phi_t(i) = phidot_t0*(time(i)-t_0) + 1/2*phiddot_max*(time(i)-t_0)^2;

        save_time1_i=i;

    elseif(time(i) >= T_1 && time(i) <= T_2)

        phiddot_t(i)=0;

        phidot_t(i)=phidot_max*0+phidot_t(save_time1_i);

        phi_t(i) = phi_t(save_time1_i) + phidot_max*(time(i)- T_1);

        save_time2_i =i;

    elseif(time(i) >= T_2 && time(i) <= T_F)

        phiddot_t(i)=-phiddot_max;

        phidot_t(i)=phidot_max*0+phidot_t(save_time2_i) - phiddot_max*(time(i) - T_2);

        phi_t(i) = phi_t(save_time2_i) + phidot_max*(time(i) - T_2) - 1/2*phiddot_max*(time(i)-T_2)^2;

        save_time3_i=i;

    else

        phiddot_t(i)=0;

        phidot_t(i)=0;

        phi_t(i)=phi_t(save_time3_i);

    end

end