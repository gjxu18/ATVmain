function [ro]=weighting_JVC
syms ro1 ro2 ro3 ro4 ro5 ro6 ro7 ro8 ro9 ro10 ro11 ro12 ro13 ro14 ro15
ro=[ro1,ro2,ro3,ro4,ro5,ro6,ro7,ro8,ro9,ro10,ro11,ro12,ro13,ro14,ro15]; 
%% 0626 roll
ro1 = 0;    %ddz
ro2 = 0;      %ddtheta
ro3 = 0 ;    %ddphi
ro4 = 11;    %zb
ro5 = 11;     %theta       
ro6 = 11;      %phi
ro7 = 1;      %tire deflection
ro8 = ro7;ro9=ro7;ro10=ro7;    
ro11 = 1;     %suspension deflection
ro12 = ro11;ro13=ro11;ro14=ro11;
ro15 = 1.5;      %u
ro = eval(ro);
%% 0626 pitch
% ro1 = 0;    %ddz
% ro2=0;      %ddtheta
% ro3 =0 ;    %ddphi
% ro4 = 11;    %zb
% ro5=11;     %theta       
% ro6=0;      %phi
% ro7=1;      %tire deflection
% ro8=ro7; ro9=ro7;ro10=ro7;    
% ro11=1;     %suspension deflection
% ro12=ro11;ro13=ro11;ro14=ro11;
% ro15=1.5;      %u
% ro=eval(ro);
%% 0625
% ro1 = 0;    %ddz
% ro2=0;      %ddtheta
% ro3 =0 ;    %ddphi
% ro4 = 1111;    %dzb
% ro5=1111;     %theta       
% ro6=1111;      %phi
% ro7=0;      %tire deflection
% ro8=ro7; ro9=ro7;ro10=ro7;    
% ro11=0;     %suspension deflection
% ro12=ro11;ro13=ro11;ro14=ro11;
% ro15=5;      %u
% ro=eval(ro);
%% lqr bump 26*26
% ro1 = 0;    %ddz
% ro2=0;      %ddtheta
% ro3 =0 ;    %ddphi
% ro4 = 1111;    %dzb
% ro5=1111;     %theta       
% ro6=1111;      %phi
% ro7=0;      %tire deflection
% ro8=ro7; ro9=ro7;ro10=ro7;    
% ro11=0;     %suspension deflection
% ro12=ro11;ro13=ro11;ro14=ro11;
% ro15=5;      %u
% ro=eval(ro);
%% 30*30 ok
% ro1 = 10^6;    %ddz
% ro2=10^5;      %ddtheta
% ro3 = 10^5;    %ddphi
% ro4 = 1;    %dzb
% ro5=1;     %theta       
% ro6=1;      %phi
% ro7=1;      %tire deflection
% ro8=ro7; ro9=ro7;ro10=ro7;    
% ro11=1;     %suspension deflection
% ro12=ro11;ro13=ro11;ro14=ro11;
% ro15=1*10^12;      %u
% ro=eval(ro);

%% dsin0617
% ro1 = 1;
% ro2=1;%tire deflection
% ro3 = 1;%suspension deflection
% ro4 = 1;%pitch&roll angle integration
% ro5=99;%controller   99     
% ro6=1;%pitch&roll angle
%% Fclass0617
% ro1 = 1;
% ro2=1;%tire deflection
% ro3 = 1;%suspension deflection
% ro4 = 1;%pitch&roll angle integration
% ro5=28;%controller       
% ro6=1;%pitch&roll angle