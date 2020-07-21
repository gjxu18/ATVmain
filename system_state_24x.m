function [A,B,D,C1,C2,C3,C34]=system_state0623(parameter)
%% ms=0；LPF；30个状态变量（u1,u1'）;jacobian；22x;y里的deg
syms mb mp ms muf mur ls1o ls2o ls3o ls4o Ca1 Ca2 Ca3 Ca4 ktf ktr...
    a b Bf Br...%Bf Br whell-track a=b=1/2wheelbase
    G Csf Csr u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 ksf ksr h Jxx Jyy Jzz Jw1 Jw2 Jw3 Jw4...
    dx1 dx2 dx3 dx4 dx5 dx6 dx7 dx8 dx9 dx10 dx11 dx12 dx13 dx14 dx15 dx16 dx17 dx18 dx19 dx20 dx21 dx22 dx23 dx24 dx25 dx26 dx27 dx28 dx29 dx30...
     x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32 x33 x34 x35 x36......
     w1 w2 w3 w4 w5 w6 w7 w8 w9 w10 w11 w12 w13 w14 u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12...
     Ip Ir hgb Fy Fx...
     fs1 fs2 fs3 fs4...%force of acuator
     ds1 ds2 ds3 ds4...
     dds1 dds2 dds3 dds4...
     dt1 dt2 dt3 dt4...
     l1 l2 l3 l4...
     dl1 dl2 dl3 dl4...
     ddl1 ddl2 ddl3 ddl4...
     mssf mssr...
     dru cfu...
     u12 u22 u32 u42...
     du12 du22 du32 du42...
%      
% mb=80;Ip=21.9;Ir=15.5;
% Bf=1.24;Br=1.24;a=0.660;b=0.660;h=0.46;
% muf=10;mur=10;
% ksf=50000;ksr=50000;Csf=2000;Csr=2000;
% ktf=200000;ktr=200000;
% msf=10;msr=10;cfu=3;dru=0.7071;

mb=parameter(1); Ip=parameter(2); Ir=parameter(3); Iz=parameter(4);
a=parameter(5); b=parameter(6); Bf=parameter(7); Br=parameter(8);
muf=parameter(9); mur=parameter(10); ktf=parameter(11); ktr=parameter(12);
ksf=parameter(13); ksr=parameter(14); Csf=parameter(15); Csr=parameter(16);
ro=parameter(17);  h=parameter(18);
msf=parameter(19); msr=parameter(20); cfu=parameter(21); dru=parameter(22);
car_speed=parameter(40);
%%
zb1=x4-a*x5+1/2*Bf*x6;
zb2=x4-a*x5-1/2*Bf*x6;
zb3=x4+b*x5+1/2*Br*x6;
zb4=x4+b*x5-1/2*Br*x6;
%%
dzb1=x1-a*x2+1/2*Bf*x3;
dzb2=x1-a*x2-1/2*Bf*x3;
dzb3=x1+b*x2+1/2*Br*x3;
dzb4=x1+b*x2-1/2*Br*x3;
%%
fs1=Csf*(x7-dzb1-x15)+ksf*(x11-zb1-x19); %u=zb-zs
fs2=Csf*(x8-dzb2-x16)+ksf*(x12-zb2-x20);
fs3=Csr*(x9-dzb3-x17)+ksf*(x13-zb3-x21);
fs4=Csr*(x10-dzb4-x18)+ksf*(x14-zb4-x22);
%% F
f1=1/mb*(fs1+fs2+fs3+fs4);%ddzb
f2=((fs3+fs4)*b-(fs1+fs2)*a)/Ip;%ddtheta
f3=((fs1-fs2)*1/2*Bf+(fs3-fs4)*1/2*Br)/Ir;%ddfai
f4=x1;f5=x2;f6=x3;
%%
f7=1/muf*(ktf*(w1-x11)-fs1);
f8=1/muf*(ktf*(w2-x12)-fs2);
f9=1/mur*(ktr*(w3-x13)-fs3);
f10=1/mur*(ktr*(w4-x14)-fs4);
f11=x7;
f12=x8;
f13=x9;
f14=x10;
f15=-2*dru*cfu*x15-cfu*cfu*x19+cfu*cfu*u1;
f16=-2*dru*cfu*x16-cfu*cfu*x20+cfu*cfu*u2;
f17=-2*dru*cfu*x17-cfu*cfu*x21+cfu*cfu*u3;
f18=-2*dru*cfu*x18-cfu*cfu*x22+cfu*cfu*u4;
f19=x15;
f20=x16;
f21=x17;
f22=x18;

%%
x=[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22];
u=[u1 u2 u3 u4];
w=[w1 w2 w3 w4];
X=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13;f14;f15;f16;f17;f18;f19;f20;f21;f22];
A=jacobian(X,x);
A = eval(A);
B=jacobian(X,u);
B = eval(B);
D=jacobian(X,w);
D = eval(D);
%% 有w的y
y1=f1;
y2=f2;
y3=f3;
y4=x4;
y5=x5/pi*180;
y6=x6/pi*180;

y7=-1*(x11-zb1-x19);
y8=-1*(x12-zb2-x20);
y9=-1*(x13-zb3-x21);
y10=-1*(x14-zb4-x22);
y11=w1-x11;
y12=w2-x12;
y13=w3-x13;
y14=w4-x14;
Y=[y1;y2;y3;y4;y5;y6;y7;y8;y9;y10;y11;y12;y13;y14];
C1=jacobian(Y,x);
C1 = eval(C1);
C2=jacobian(Y,u);
C2 = eval(C2);
C3=jacobian(Y,w);
C3 = eval(C3);
%% w在里的y
x=[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26];
y11=x23-x11;
y12=x24-x12;
y13=x25-x13;
y14=x26-x14;
Y=[y1;y2;y3;y4;y5;y6;y7;y8;y9;y10;y11;y12;y13;y14];
C34=jacobian(Y,x);
C34 = eval(C34);



% u12=zb1-x11;
% u22=zb2-x12;
% u32=zb3-x13;
% u42=zb4-x14;
% zr1=zb1-x11-x19-x27;
% zr2=zb2-x12-x20-x28;
% zr3=zb3-x13-x21-x29;
% zr4=zb4-x14-x22-x30;

