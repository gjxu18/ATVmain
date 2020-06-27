function [A,B,D,C1,C2,C3,C34]=system_state_30x(parameter)
%% ms=0；LPF；30个状态变量（u1,u1'）;jacobian；
syms Fxs1 Fxs2 Fxs3 Fxs4 Fys1 Fys2 Fys3 Fys4 Fzs1 Fzs2 Fzs3 Fzs4... 
    Fxg1 Fxg2 Fxg3 Fxg4 Fyg1 Fyg2 Fyg3 Fyg4 Fzg1 Fzg2 Fzg3 Fzg4...
    mb mp ms muf mur ls1o ls2o ls3o ls4o Ca1 Ca2 Ca3 Ca4 ktf ktr...
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
ma=10;
mb=parameter(1); Ip=parameter(2); Ir=parameter(3); Iz=parameter(4);
a=parameter(5); b=parameter(6); Bf=parameter(7); Br=parameter(8);
muf=parameter(9); mur=parameter(10); ktf=parameter(11); ktr=parameter(12);
ksf=parameter(13); ksr=parameter(14); Csf=parameter(15); Csr=parameter(16);
ro=parameter(17);  h=parameter(18);
msf=parameter(19); msr=parameter(20); cfu=parameter(21); dru=parameter(22);
car_speed=parameter(40);
%%
% ddx23=-2*dru*cfu*dx23-cfu*cfu*x23+cfu*cfu*x27;
% ddx24=-2*dru*cfu*dx24-cfu*cfu*x24+cfu*cfu*x28;
% ddx25=-2*dru*cfu*dx25-cfu*cfu*x25+cfu*cfu*x29;
% ddx26=-2*dru*cfu*dx26-cfu*cfu*x26+cfu*cfu*x30;
%%
zb1=x4-a*x5+1/2*Bf*x6;
zb2=x4-a*x5-1/2*Br*x6;
zb3=x4+b*x5+1/2*Br*x6;
zb4=x4+b*x5-1/2*Br*x6;
%%
dzb1=x1-a*x2+1/2*Bf*x3;
dzb2=x1-a*x2-1/2*Bf*x3;
dzb3=x1+b*x2+1/2*Br*x3;
dzb4=x1+b*x2-1/2*Br*x3;
%%
u12=zb1-x11;
u22=zb2-x12;
u32=zb3-x13;
u42=zb4-x14;
%%
du12=dzb1-x7;
du22=dzb2-x8;
du32=dzb3-x9;
du42=dzb4-x10;
%%
ddu12=cfu*cfu*x27-2*dru*cfu*du12-cfu*cfu*u12;
ddu22=cfu*cfu*x28-2*dru*cfu*du22-cfu*cfu*u22;
ddu32=cfu*cfu*x29-2*dru*cfu*du32-cfu*cfu*u32;
ddu42=cfu*cfu*x30-2*dru*cfu*du42-cfu*cfu*u42;
%%
% zb1=zs1+zw1+zr1+l1;
% zb2=zs2+zw2+zr2+l2;
% zb3=zs3+zw3+zr3+l3;
% zb4=zs4+zw4+zr4+l4;

%%
fs1=Csf*(x15-x7)+ksf*(x19-x11);
fs2=Csf*(x16-x8)+ksf*(x20-x12);
fs3=Csr*(x17-x9)+ksf*(x21-x13);
fs4=Csr*(x18-x10)+ksf*(x22-x14);
%% 
% Fx=w5;Fy=w6;
%% 
External_Fxy=[Fx Fy];
%% F
f1=1/mb*(fs1+fs2+fs3+fs4);%ddzb
f2=((fs3+fs4)*b-(fs1+fs2)*a+Fx.*h)/Ip;%ddtheta
f3=((fs1-fs2)*1/2*Bf+(fs3-fs4)*1/2*Br+Fy.*h)/Ir;%ddfai
f4=x1;f5=x2;f6=x3;
%%
dx1=f1;dx2=f2;dx3=f3;
ddzb1=dx1-a*dx2+1/2*Bf*dx3;
ddzb2=dx1-a*dx2-1/2*Bf*dx3;
ddzb3=dx1+b*dx2+1/2*Br*dx3;
ddzb4=dx1+b*dx2-1/2*Br*dx3;
%%
ddzs1=1/msf*(fs1-ma*ddu12);
ddzs2=1/msf*(fs2-ma*ddu22);
ddzs3=1/msr*(fs3-ma*ddu32);
ddzs4=1/msr*(fs4-ma*ddu42);
%%
f7=ddzs1;
f8=ddzs2;
f9=ddzs3;
f10=ddzs4;
f11=x7;
f12=x8;
f13=x9;
f14=x10;
f15=1/muf*(Csf*(x7-x15)+ksf*(x11-x19)+ktf*(w1-x19));%可以加上u''+模型验证
f16=1/muf*(Csf*(x8-x16)+ksf*(x12-x20)+ktf*(w2-x20));
f17=1/mur*(Csr*(x9-x17)+ksr*(x13-x21)+ktr*(w3-x21));
f18=1/mur*(Csr*(x10-x18)+ksr*(x14-x22)+ktr*(w4-x22));
f19=x15;
f20=x16;
f21=x17;
f22=x18;
f23=-2*dru*cfu*x23-cfu*cfu*x27+cfu*cfu*u1;
f24=-2*dru*cfu*x24-cfu*cfu*x28+cfu*cfu*u2;
f25=-2*dru*cfu*x25-cfu*cfu*x29+cfu*cfu*u3;
f26=-2*dru*cfu*x26-cfu*cfu*x30+cfu*cfu*u4;
f27=x23;
f28=x24;
f29=x25;
f30=x26;
%%
x=[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30];
u=[u1 u2 u3 u4];
w=[w1 w2 w3 w4];
X=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13;f14;f15;f16;f17;f18;f19;f20;f21;f22;f23;f24;f25;f26;f27;f28;f29;f30];
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
y5=x5;
y6=x6;
y7=w1-x19;
y8=w2-x20;
y9=w3-x21;
y10=w4-x22;
y11=x19-zb1;
y12=x20-zb2;
y13=x21-zb3;
y14=x22-zb4;
Y=[y1;y2;y3;y4;y5;y6;y7;y8;y9;y10;y11;y12;y13;y14];
C1=jacobian(Y,x);
C1 = eval(C1);
C2=jacobian(Y,u);
C2 = eval(C2);
C3=jacobian(Y,w);
C3 = eval(C3);
%% w在里的y
x=[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32 x33 x34];
y7=x31-x19;
y8=x32-x20;
y9=x33-x21;
y10=x34-x22;
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

