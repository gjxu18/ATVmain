function [rho_kal]=weighting_Kalman
syms ro1 ro2 ro3 ro4 ro5 ro6 ro7 ro8 ro9 ro10 ro11 ro12 ro13 ro14 ro15 ro16 ro17 ro18 ro19 ro20
ro=[ro1,ro2,ro3,ro4,ro5,ro6,ro7,ro8,ro9,ro10,ro11,ro12,ro13,ro14,ro15,ro16,ro17,ro18,ro19,ro20]; 
%%  F clss road 0719
ro1 = 9;    %Q1
ro2 = 0;      %Q2
ro3 = 9 ;    %Q3
ro4 = 0;    %Q4
ro5 = 0.000001;     %R       
ro6 = 111;      %phi
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
rho_kal = eval(ro);
