function [rho]=weighting_MPC
syms ro1 ro2 ro3 ro4 ro5 ro6 ro7 ro8 ro9 ro10 ro11 ro12 ro13 ro14 ro15 ro16 ro17 ro18 ro19 ro20
ro=[ro1,ro2,ro3,ro4,ro5,ro6,ro7,ro8,ro9,ro10,ro11,ro12,ro13,ro14,ro15,ro16,ro17,ro18,ro19,ro20]; 
%% 0708 mpc bump
ro1 = 0;    %ddz
ro2 = 0;      %ddtheta
ro3 = 0 ;    %ddphi
ro4 = 1;    %zb
ro5 = 1;     %theta       
ro6 = 1;      %phi
ro7 = 0;      %tire deflection
ro8 = ro7;ro9=ro7;ro10=ro7;    
ro11 = 0;     %suspension deflection
ro12 = ro11;ro13=ro11;ro14=ro11;
ro15 = 1;      %u1
ro16 = 1;
ro17 = 1;
ro18 = 1;       %u4

ro19 = 25;      %p H
ro20 = 5;       %m H
rho = eval(ro);