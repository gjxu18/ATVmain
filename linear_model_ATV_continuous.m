% % 20200627 整理：包含共6个函数：参数，状态空间2个，权重，计算QR，lqr
clc;
clear;
syms x
%% vehicle paremeter
parameter=atv_parameters;
mb=parameter(1); Ip=parameter(2); Ir=parameter(3); Iz=parameter(4);
a=parameter(5); b=parameter(6); Bf=parameter(7); Br=parameter(8);
muf=parameter(9); mur=parameter(10); ktf=parameter(11); ktr=parameter(12);
ksf=parameter(13); ksr=parameter(14); Csf=parameter(15); Csr=parameter(16);
ro=parameter(17);  h=parameter(18);
msf=parameter(19); msr=parameter(20); cfu=parameter(21); dru=parameter(22);
car_speed=parameter(40);
%% sim set
Tstop = 15;            %stoptime
Ts = 0.01;           %1/(10*car_speed); %sample time
tout = 0:Ts:Tstop;     %time
xstop = Tstop/Ts+1;    %x' number    
timedelay=(a+b)/car_speed;  %front-rear time
xdelay=fix(timedelay/Ts); %取整   %front-rear number
%% dotx=Ax+Bu+Dw
[Ap,Bp,Dp,Cp,C2p,C3p,C26]=system_state_24x(parameter);
% [Ap,Bp,Dp,Cp,C2p,C3p,C26]=system_state_30x(parameter);
[SXp,SYp]=size(Ap);
%% controllability test
disp('controllability test')
conte=[];
for i=1:SXp
    conte=[conte Ap^(i-1)*Bp];
end
rank(conte)
if rank(conte) == SXp
disp('The system is controllable.')
else 
disp('The system is not controllable.')
end

%% augmented matrix
% S=zeros(2,SXp);
% S(1,5)=1;S(2,6)=1;
% As=[Ap zeros(30,2);
%     S zeros(2,2)];
% Bs=[Bp;zeros(2,4)];
% Ds=[Dp;zeros(2,4)];
% C2=zeros(14,2);
% Cs=[Cp C2];
%% 离散化
% Asd=expm(As.*Ts);
% fun=@(x)expm(As.*x);
% Bsd=integral(fun,0,Ts,'ArrayValued',true)*Bs;
% Dsd=integral(fun,0,Ts,'ArrayValued',true)*Ds;
%% 等级路面 dotx=Aw*x+Iw*w
av=2*pi*0.01*car_speed;n0=0.1;
G0=16284*10^(-6); %F级路面16284 D级1024 C级256 E级4096
Aw=zeros(4,4);Aw(1,1)=-av;Aw(2,2)=-av;Aw(3,3)=-av;Aw(4,4)=-av;
Iw=eye(4);Iw(1,1)=2*pi*n0*sqrt(G0*car_speed);Iw(2,2)=2*pi*n0*sqrt(G0*car_speed);
Iw(3,3)=2*pi*n0*sqrt(G0*car_speed);Iw(4,4)=2*pi*n0*sqrt(G0*car_speed);
Iw_0=Iw^(-1);
noise1=wgn(xstop,1,30,'linear');
y_n=var(noise1);
noise2=[zeros(xdelay,1);noise1];

noise_f=[tout' noise1];
noise_r=[tout' noise2(1:xstop)];

figure('name','class road')
subplot(2,2,1)
plot(tout,noise_f(:,2));
subplot(2,2,2)
plot(tout,noise_f(:,2));
subplot(2,2,3)
plot(tout,noise_r(:,2));
subplot(2,2,4)
plot(tout,noise_r(:,2));
%% sin波形路面
a_max=0.1342;
sin_f=a_max*sin(2*pi/6*tout');
sin_r=[zeros(xdelay,1);sin_f];
sin_0=0.*tout';

sin_f=[tout' sin_f];
sin_r=[tout' sin_r(1:xstop)];
sin_0=[tout' sin_0];

figure('name','sin road')
subplot(2,2,1)
plot(tout,sin_f(:,2));
subplot(2,2,2)
plot(tout,sin_f(:,2));
subplot(2,2,3)
plot(tout,sin_r(:,2));
subplot(2,2,4)
plot(tout,sin_r(:,2));

%% Bump Road
Abump=0.1342; %坑、包的幅值
Lbump=3;    %坑、包的长度
t_hc=3;     %缓冲时间
xbump_t=0:Ts:(Lbump/car_speed);
xbump_size=size(xbump_t,2);
xbump=Abump/2*(1-cos(2*pi*car_speed/Lbump*xbump_t));
xxx=-1.*xbump;
xbump_fl=[zeros(1,3/Ts) xbump -1.*xbump zeros(1,xstop-2*xbump_size-3/Ts)];
xbump_rl=[zeros(1,3/Ts+xdelay) xbump -1.*xbump zeros(1,xstop-2*xbump_size-3/Ts-xdelay)];
xbump_fr=zeros(1,xstop);
xbump_rr=zeros(1,xstop);

figure('name','bump road')
subplot(2,2,1)
plot(tout,xbump_fl);
subplot(2,2,2)
plot(tout,xbump_fr);
subplot(2,2,3)
plot(tout,xbump_rl);
subplot(2,2,4)
plot(tout,xbump_rr);

bump_fl=[tout' xbump_fl'];
bump_rl=[tout' xbump_rl'];
bump_fr=[tout' xbump_fr'];
bump_rr=[tout' xbump_rr'];

bump4=[xbump_fl;xbump_fr;xbump_rl;xbump_rr];
%% 3D road y以mm为单位
% ybump=0:1:2000;
% zbump=zeros(900+1,xbump_size);
% for i=ybump  
%     if i<300
%         zbump(i+1,:)=(Abump*i/300)/2*(1-cos(2*pi*car_speed/Lbump*xbump_t));
%     elseif i>=300 && i<600
%         zbump(i+1,:)=(Abump)/2*(1-cos(2*pi*car_speed/Lbump*xbump_t));
%     elseif i>=600 && i<900 
%         zbump(i+1,:)=(900-i)/300*(Abump)/2*(1-cos(2*pi*car_speed/Lbump*xbump_t));
%     elseif i>=900 && i<1100
%         zbump(i+1,:)=0.*xbump_t;
%     elseif i>=1100 && i<1400
%         zbump(i+1,:)=(Abump*(i-1100)/300)/2*(1-cos(2*pi*car_speed/Lbump*xbump_t));
%     elseif i>=1400 && i<1700
%         zbump(i+1,:)=(Abump)/2*(1-cos(2*pi*car_speed/Lbump*xbump_t));
%     else 
%         zbump(i+1,:)=(2000-i)/300*(Abump)/2*(1-cos(2*pi*car_speed/Lbump*xbump_t));
%     end
% end
% zbump=[zeros(2001,3/Ts) zbump -1.*zbump zeros(2001,xstop-2*xbump_size-3/Ts)];
% zbump=[zeros(300,1501);zbump;zeros(300,1501)];
% ybump1=0:1:2600;
% figure ('name','3d road')
% surf(tout,ybump1,zbump);
% mesh(tout,ybump1,zbump);
% zbump_origin=zbump(:,1:1201);
%% dotx=Ax+B(u w) 可行很棒！未验证
% A=[As zeros(32,4);zeros(4,36)];
% B=[Bs Ds;
%     zeros(4,8)];
% C=[Cp zeros(14,6)];
%% +4
A=[Ap Dp;
    zeros(4,SXp) Aw];
B=[Bp;zeros(4,4)];
D=[zeros(SXp,4);Iw];
C=C26;
[SX,SY]=size(A);
%% control matrix
[ro]=weighting_JVC;
%% lqr
% [Q,R]=optimal_control_matrix_JVC_no_w(parameter,A,B,SX,SY,ro);
Q=eye(14);R=eye(4);
for i=1:14
    Q(i,i)=ro(i);
end
for i=1:4
    R(i,i)=ro(15);
end
[Abar,Bbar,Cbar,T,k]=ctrbf(A,B,C);
% [P,P12,P11]=transfer_uncontrollable_matrix_continuous_time(A,B,C,D,k,R,Q);
Q=C'*Q*C;
[K,S,P]=lqr(A,B,Q,R);
% K=R^(-1)*B'*P;
Ac=A-B*K;
% Pk=[P11 P12];      %论文和上边一样的
% K1=R^(-1)*Bp'*Pk;
% Ac=A-B*K;
cun=size(Ac);
Cu=zeros(4,cun(1)); %得到u'
Cu(1,19)=1;
Cu(2,20)=1;
Cu(3,21)=1;
Cu(4,22)=1;
%% lqg
Dn=zeros(14,8);
% Qn=35.*eye(4);   0704
% Rn=0.001.*eye(14);
Qn=500.*eye(4); 
Rn=0.0001.*eye(14);
Nn=zeros(4,14);
sys=ss(A,[B D],C,Dn);
[kest,L,P]=kalman(sys,Qn,Rn,Nn);
regulator=lqgreg(kest,K);
%% mpc
Cmpc=zeros(3,26);
Cmpc(1,4)=1; %其它有用
Cmpc(2,5)=1;
Cmpc(3,6)=1; 
% Cmpc(4,4)=1;
% Cmpc(5,5)=1;
% Cmpc(6,6)=1;
% Cmpc(1,4)=1;
% Cmpc(2,5)=1;
% Cmpc(3,6)=1;
refer=zeros(3,1);
Cmpcu=zeros(4,22); %得到u'
Cmpcu(1,19)=1;
Cmpcu(2,20)=1;
Cmpcu(3,21)=1;
Cmpcu(4,22)=1;
%% kalman
Dn=zeros(14,8);
% Qn=35.*eye(4);   0704
% Rn=0.001.*eye(14);
Qn=1.*eye(4); 
Rn=0.00001.*eye(14);
Nn=zeros(4,14);
sys=ss(A,[B D],C,Dn);
[kest,L,P]=kalman(sys,Qn,Rn,Nn);

%%
% Q=eye()
ob=obsv(A,C);
obb=rank(ob);
[Abar,Bbar,Cbar,T,k] = obsvf(A,B,C);
sum(k)
Ckalman=zeros(9,14);
Ckalman(1,1)=1;
Ckalman(2,2)=1;
Ckalman(3,3)=1;
Ckalman(4,5)=1;
Ckalman(5,6)=1;
Ckalman(6,7)=1;
Ckalman(7,8)=1;
Ckalman(8,9)=1;
Ckalman(9,10)=1;
Ckalman=Ckalman*C;
Ckalman4=zeros(8,26);
Ckalman4(1,15)=1;
Ckalman4(2,16)=1;
Ckalman4(3,17)=1;
Ckalman4(4,18)=1;
Ckalman4(5,19)=1;
Ckalman4(6,20)=1;
Ckalman4(7,21)=1;
Ckalman4(8,22)=1;
Ckalman=[Ckalman;Ckalman4];
%% MPC bump
%******连续模型*******%
Ac=Ap;Bcu=Bp;Bcd=Dp;
Cmb=zeros(4,14);
Cmb(1,1)=1;
Cmb(2,4)=1;
Cmb(3,5)=1;
Cmb(4,6)=1;
Cmb=Cmb*Cp;
Cc=Cmb;                         %输出矩阵
Cm=Cmb;                         %测量矩阵（两者可能不同）
%******离散化*******%
Ad=expm(Ac*Ts);
fun=@(x)expm(Ac*x);
Bu=integral(fun,0,Ts,'ArrayValued',true)*Bcu;
Bd=integral(fun,0,Ts,'ArrayValued',true)*Bcd;

% [m,n] = size(Ac); %#ok<ASGLU>
% [m,nb] = size(Bcd); %#ok<ASGLU>
% s = expm([[Ac Bcd]*Ts; zeros(nb,n+nb)]);
% Ad = s(1:n,1:n);
% Bd = s(1:n,n+1:n+nb);
%*********增量状态空间模型********%
[xm,xm]=size(Ad);[xm,um]=size(Bu);[xm,dm]=size(Bd);
[ym,xm]=size(Cc);
%**********初始化*************%
[rho]=weighting_MPC;
p=rho(16);m=rho(17);
rho_y=[rho(1) rho(4:6)];rho_u=rho(15).*ones(1,4);
refer=zeros(ym,1);                      %reference自定义  修改
Refer=zeros(ym*p,xstop);
[Kmpc,Sx,I,Sd]=Model_Predictive_Control(Ad,Bu,Bd,Cc,Ts,p,m,rho_y,rho_u);
%**********求delta d road**********%
road4=bump4;
dk_road=zeros(4,xstop);
for i=2:xstop
    dk_road(:,i)=road4(:,i)-road4(:,i-1);
end
%*************求测量值ym***********%
x=zeros(xm,xstop);
yc=zeros(ym,xstop);
Ym=zeros(ym,xstop);
% x(:,i)=Cm/Ym(:,i);
%*********计算误差与控制量变化量***********%
Ep=zeros(ym*p,xstop);
dx=zeros(xm,xstop);
du=zeros(um,xstop);
for i=2:xstop-1
%     Ep(:,1+1)=Refer(:,1+1)-Sx*dx(:,1)-I*yc(:,1)-Sd*dk_road(:,1);
%     du(:,1)=Kmpc*Ep(:,1+1);
%     dx(:,1+1)=Ad*dx(:,1)+Bu*du(:,1)+Bd*dk_road(:,1);
%     yc(:,1+1)=Cc*dx(:,1)+zeros(ym,1);
    
    Ep(:,i+1)=Refer(:,i+1)-Sx*dx(:,i)-I*yc(:,i)-Sd*dk_road(:,i);
    du(:,i)=Kmpc*Ep(:,i+1);
    dx(:,i+1)=Ad*dx(:,i)+Bu*du(:,i)+Bd*dk_road(:,i);
    yc(:,i+1)=Cc*dx(:,i)+yc(:,i-1);
end
for i=1:xstop-1
    x(:,i+1)=x(:,i)+dx(:,i+1);
end
%***********被动比较***********%
xp=zeros(xm,xstop);
yp=zeros(ym,xstop);
for i=1:xstop
    xp(:,i+1)=Ad*xp(:,i)+Bd*road4(:,i);
    yp(:,i)=Cc*xp(:,i);
end

figure('name','compare')
subplot(2,2,1)
plot(tout,yp(1,:),tout,yc(1,:));
subplot(2,2,2)
plot(tout,yp(2,:),tout,yc(2,:));
subplot(2,2,3)
plot(tout,yp(3,:),tout,yc(3,:));
subplot(2,2,4)
plot(tout,yp(4,:),tout,yc(4,:));

figure('name','u_mpc')
subplot(2,2,1)
plot(tout,x(19,:));
subplot(2,2,2)
plot(tout,x(20,:));
subplot(2,2,3)
plot(tout,x(21,:));
subplot(2,2,4)
plot(tout,x(22,:));