function paramter=atv_parameters
        Carf=33*1000;Calf=33*1000;Calr=33*1000;Carr=33*1000;
%         Carf=51.285*1000;Calf=51.285*1000;Calr=51.285*1000;Carr=51.285*1000;%ПаµИ
%         Carf=11*1000;Calf=11*1000;Calr=11*1000;Carr=11*1000;        
        Ctrf=35000;Ctlf=35000;Ctlr=35000;Ctrr=35000;
        Jtrf=2.1;Jtlf=2.1;Jtlr=2.1;Jtrr=2.1;
        hrcf=0.5;hrcr=0.5;Cf=1.5;Cr=1.5;
        
  G=9.8; car_speed=1;
mb=142.8  %80-110 adams34.2   dianchi50kg  chejia30.36kg
Ip=8.145;  %10-21.9
Ir=9.55;Iz=8;  %8-15.5
Bf=1.3;Br=1.3;  %0.8-1.24
a=0.6099;b=0.6099;
h=0.46;ro=0.285;
muf=22;mur=22;
ksf=17000;ksr=17000;%18000  
Csf=1000;Csr=1000;  %548
ktf=100000;ktr=100000;%200000      175500
msf=6.6;msr=6.6;
cfu=30*pi;dru=0.7071;%0.7071;%cut-off-frequency damping-ratio
%% 0626      
% G=9.8; car_speed=1;
% mb=142.8  %80-110 adams34.2   dianchi50kg  chejia30.36kg
% Ip=8.145;  %10-21.9
% Ir=9.55;Iz=8;  %8-15.5
% Bf=1.3;Br=1.3;  %0.8-1.24
% a=0.6099;b=0.6099;
% h=0.46;ro=0.285;
% muf=22;mur=22;
% ksf=170000;ksr=170000;%18000  
% Csf=1000;Csr=1000;  %548
% ktf=100000;ktr=100000;%200000      175500
% msf=6.6;msr=6.6;
% cfu=2*pi;dru=Csf/2/(sqrt(ksf*mb/4));%0.7071;%cut-off-frequency damping-ratio
%         
%%  
paramter=[mb Ip Ir Iz...
            a    b   Bf   Br...
            muf mur ktf ktr...
            ksf ksr Csf Csr...
            ro  h...
            msf msr cfu dru...
            Carf Calf Calr Carr...
            Ctrf Ctlf Ctlr Ctrr...
            Jtrf Jtlf Jtlr Jtrr...
            hrcf hrcr Cf   Cr...
            G car_speed];