%% 20200718 离散 kalman状态估计 ok！ 
%******连续模型*******%
Ac=A;Bcu=B;Bcd=D;
Cmb=zeros(3,14);
Cmb(1,1)=1;
Cmb(2,5)=1;
Cmb(3,6)=1;
Cmb(4,6)=1;
Cmb=Cmb*C;
Cc=Cmb;                         %输出矩阵
Cm=Cmb;                         %测量矩阵（两者可能不同）
%******离散化*******%
Ad=expm(Ac*Ts);
fun=@(x)expm(Ac*x);
Bu=integral(fun,0,Ts,'ArrayValued',true)*Bcu;
Bd=integral(fun,0,Ts,'ArrayValued',true)*Bcd;
%*********状态空间模型********%
[xm,xm]=size(Ad);[xm,um]=size(Bu);[xm,dm]=size(Bd);
[ym,xm]=size(Cc);Du=zeros(ym,um+dm);
%**********初始化*************%
[rho_kal]=weighting_Kalman;
Qk=zeros(dm,dm);
for i=1:dm
    Qk(i,i)=rho_kal(i);
end
Rk=0.00001.*eye(ym);Nn=zeros(dm,ym);
%**********求delta d road**********%
road4=bump4;
% road4=sin4;
w4=zeros(dm,xstop);
for i=1:xstop
    if i==xstop
        w4(:,i)=Bd(xm-4+1:xm,1:dm)^(-1)*(road4(:,i)-Ad(xm-4+1:xm,xm-4+1:xm)*road4(:,i));
    else
        w4(:,i)=Bd(xm-4+1:xm,1:dm)^(-1)*(road4(:,i+1)-Ad(xm-4+1:xm,xm-4+1:xm)*road4(:,i));
    end
end
% w4=noise4;
w42sim=[tout' w4'];
dk_w=zeros(4,xstop);
for i=2:xstop
    dk_w(:,i)=w4(:,i)-w4(:,i-1);
end
%*************求测量值ym***********%
x=zeros(xm,xstop);
yc=zeros(ym,xstop);
Ym=zeros(ym,xstop);
noise_vk=wgn(ym,xstop,0.000001,'linear');
u=zeros(um,xstop);
P_k=0.01.*eye(xm);
x_hat=zeros(xm,xstop);
xp=zeros(xm,xstop);
yp=zeros(ym,xstop);
for i=1:xstop-1
    %***********被动***********%
    xp(:,i+1)=Ad*xp(:,i)+Bd*w4(:,i)+Bu*u(:,i);
    yp(:,i)=Cc*xp(:,i);
    Ym(:,i)=yp(:,i)+noise_vk(:,i);
    %**************Kalman_***********%        
        K_kalman=P_k*Cm'/(Cm*P_k*Cm'+Rk);
        x_hat(:,i)=x(:,i)+K_kalman*(Ym(:,i)-Cm*x(:,i));
        Pk=(eye(xm)-K_kalman*Cm)*P_k;
        P_k=Ad*Pk*Ad'+Bd*Qk*Bd';
        x(:,i+1)=Ad*x_hat(:,i)+Bu*u(:,i);%+Bd*w4(:,i); 
        yc(:,i)=Cc*x_hat(:,i);  
end
Croad=zeros(dm,xm);
Croad(1,23)=1;
Croad(2,24)=1;
Croad(3,25)=1;
Croad(4,26)=1;
figure('name','road compare')
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
plot(tout,Ym(1,:),tout,yc(1,:));
title('\alpha_w');
subplot(2,2,2)
plot(tout,Ym(2,:),tout,yc(2,:));
title('z');
subplot(2,2,3)
plot(tout,Ym(3,:),tout,yc(3,:));
title('\theta');
subplot(2,2,4)
% plot(tout,Ym(4,:),tout,yc(4,:));
% title('\phi');
