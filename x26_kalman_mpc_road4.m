%% MPC x26 无约束 有估计 估计路面
%******连续模型*******%
Ac=A;Bcu=B;Bcd=D;
Cmb=zeros(4,14);
Cmb(1,1)=1;
Cmb(2,4)=1;
Cmb(3,5)=1;
Cmb(4,6)=1;
Cmb=Cmb*C;
Cc=Cmb;                         %输出矩阵
Cm=Cmb;                         %测量矩阵（两者可能不同）
%******离散化*******%
Ad=expm(Ac*Ts);
fun=@(x)expm(Ac*x);
Bu=integral(fun,0,Ts,'ArrayValued',true)*Bcu;
Bd=integral(fun,0,Ts,'ArrayValued',true)*Bcd;
%*********增量状态空间模型********%
[xm,xm]=size(Ad);[xm,um]=size(Bu);[xm,dm]=size(Bd);
[ym,xm]=size(Cc);Du=zeros(ym,um+dm);
%**********初始化*************%
[rho]=weighting_MPC;
p=rho(16);m=rho(17);
rho_y=[rho(1) rho(4:6)];rho_u=rho(15).*ones(1,4);
refer=zeros(ym,1);                      %reference自定义  修改
Refer=zeros(ym*p,xstop);
[Kmpc,Sx,I,Sd]=Model_Predictive_Control(Ad,Bu,Bd,Cc,Ts,p,m,rho_y,rho_u);
[rho_kal]=weighting_Kalman;
Qk=zeros(dm,dm);
for i=1:4
    Qk(i,i)=rho_kal(i);
end
Rk=0.00001.*eye(ym);Nn=zeros(dm,ym);
%**********求delta d road w4**********%
road4=bump4;
for i=1:xstop
    if i==xstop
        w4(:,i)=Bd(xm-4+1:xm,1:dm)^(-1)*(road4(:,i)-Ad(xm-4+1:xm,xm-4+1:xm)*road4(:,i));
    else
        w4(:,i)=Bd(xm-4+1:xm,1:dm)^(-1)*(road4(:,i+1)-Ad(xm-4+1:xm,xm-4+1:xm)*road4(:,i));
    end
end
% w4=noise4;
dk_w=zeros(4,xstop);
for i=2:xstop
    dk_w(:,i)=w4(:,i)-w4(:,i-1);
end
%*************求测量值ym***********%
x=zeros(xm,xstop);
yc=zeros(ym,xstop);
Xm=zeros(xm,xstop);
Ym=zeros(ym,xstop);
noise_vk=wgn(4,xstop,0.0000001,'linear');
%*********计算误差与控制量变化量***********%
Ep=zeros(ym*p,xstop);
dx=zeros(xm,xstop);
du=zeros(um,xstop);
u=zeros(um,xstop);
P_k=0.01.*eye(xm);
x_hat=zeros(xm,xstop);
dx_hat=zeros(xm,xstop);
yc_hat=zeros(ym,xstop);
xp=zeros(xm,xstop);
yp=zeros(ym,xstop);
for i=1:xstop-1
    %***********得到测量值***********%
    Ym(:,i)=Cc*Xm(:,i)+noise_vk(:,i);
    %**************Kalman_***********%   
    K_kalman=P_k*Cm'/(Cm*P_k*Cm'+Rk);
    x_hat(:,i)=x(:,i)+K_kalman*(Ym(:,i)-Cm*x(:,i));
    Pk=(eye(xm)-K_kalman*Cm)*P_k;
    P_k=Ad*Pk*Ad'+Bd*Qk*Bd';
    x(:,i+1)=Ad*x_hat(:,i)+Bu*u(:,i);%+Bd*w4(:,i); 
    yc(:,i)=Cc*x_hat(:,i);  
    if i==1
        dx_hat(:,i)=x_hat(:,i)-x_hat(:,i);
    else
        dx_hat(:,i)=x_hat(:,i)-x_hat(:,i-1);
    end
    %*************计算误差***********%
    Ep(:,i+1)=Refer(:,i+1)-Sx*dx_hat(:,i)-I*yc(:,i)-Sd*dk_road(:,i);
    du(:,i)=Kmpc*Ep(:,i+1);
    if i==1
        u(:,i)=u(:,i)+du(:,i);%u(-1)=0
    else
        u(:,i)=u(:,i-1)+du(:,i);
    end 
    Xm(:,i+1)=Ad*x_hat(:,i)+Bd*w4(:,i)+Bu*u(:,i);
        %***********被动比较***********%
    xp(:,i+1)=Ad*xp(:,i)+Bd*w4(:,i);%+Bu*u(:,i);
    yp(:,i)=Cc*xp(:,i);
end


figure('name','road')
subplot(2,2,1)
plot(tout,xp(23,:),tout,x_hat(23,:));
subplot(2,2,2)
plot(tout,xp(24,:),tout,x_hat(24,:));
subplot(2,2,3)
plot(tout,xp(25,:),tout,x_hat(25,:));
subplot(2,2,4)
plot(tout,xp(26,:),tout,x_hat(26,:));

figure('name','compare')
subplot(2,2,1)
plot(tout,yp(1,:),tout,yc(1,:));
title('\alpha_w');
subplot(2,2,2)
plot(tout,yp(2,:),tout,yc(2,:));
title('z');
subplot(2,2,3)
plot(tout,yp(3,:),tout,yc(3,:));
title('\theta');
subplot(2,2,4)
plot(tout,yp(4,:),tout,yc(4,:));
title('\phi');

figure('name','u_mpc')
subplot(2,2,1)
plot(tout,x_hat(19,:));
subplot(2,2,2)
plot(tout,x_hat(20,:));
subplot(2,2,3)
plot(tout,x_hat(21,:));
subplot(2,2,4)
plot(tout,x_hat(22,:));