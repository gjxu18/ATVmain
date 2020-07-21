%% MPC x22 有约束 无估计
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
Cb=Cmb;
%******离散化*******%
Ad=expm(Ac*Ts);
fun=@(x)expm(Ac*x);
Bu=integral(fun,0,Ts,'ArrayValued',true)*Bcu;
Bd=integral(fun,0,Ts,'ArrayValued',true)*Bcd;
%*********增量状态空间模型********%
[xm,xm]=size(Ad);[xm,um]=size(Bu);[xm,dm]=size(Bd);
[ym,xm]=size(Cc);[bm,xm]=size(Cb);
%**********初始化*************%
[rho]=weighting_MPC;
p=rho(19);m=rho(20);
rho_y=[rho(1) rho(4:6)];rho_u=[rho(15) rho(16) rho(17) rho(18)] ;
refer=zeros(ym,1);                      %reference自定义  修改
Refer=zeros(ym*p,xstop);
[Kg,Sx,I,Sd,Sxb,Sub,Sdb,Ib,H,Cu]=MPC_constraint(Ad,Bu,Bd,Cc,Cb,Ts,p,m,rho_y,rho_u);
%**********求delta d road**********%
road4=bump4;
% road4=sin4;
dk_road=zeros(4,xstop);
for i=2:xstop
    dk_road(:,i)=road4(:,i)-road4(:,i-1);
end
%*************求测量值ym***********%
x=zeros(xm,xstop);
yc=zeros(ym,xstop);
%*********计算误差与控制量变化量***********%
Ep=zeros(ym*p,xstop);
G=zeros(um*m,xstop);
dx=zeros(xm,xstop);
du=zeros(um,xstop);
u=zeros(um,xstop);
xp=zeros(xm,xstop);
yp=zeros(ym,xstop);
%*************约束*****************%
bk=zeros(um*m*2+um*m*2+bm*p*2,xstop);
bk=zeros(um*m*2+um*m*2,xstop);
for i=2:xstop-1 
    if i==1
        dx(:,i)=x(:,i)-x(:,i);%x(-1)=0
    else
        dx(:,i)=x(:,i)-x(:,i-1);
    end
    yc(:,i)=Cc*x(:,i);
    
    Ep(:,i+1)=Refer(:,i+1)-Sx*dx(:,i)-I*yc(:,i)-Sd*dk_road(:,i);
    G(:,i+1)=Kg*Ep(:,i+1);
    for j=1:um*m*2+um*m*2+bm*p*2
        if j<=m*um
            bk(j,i+1)=-rho(21);
        elseif j<=2*m*um && j>m*um
            bk(j,i+1)=rho(21);
        elseif j>= 2*m*um+1 && j<=2*m*um+m
            bk(2*m*um+4*(j-2*m*um)-3:2*m*um+4*(j-2*m*um),i+1)=u(:,i-1)-rho(22).*ones(um,1);
        elseif j>= 3*um*m+1 && j<=3*m*um+m
            bk(3*m*um+4*(j-3*m*um)-3:3*m*um+4*(j-3*m*um),i+1)=rho(23).*ones(um,1)-u(:,i-1);
%         elseif j<=4*um*m+bm*p && j>4*um*m      %无Y约束
%             bk(j,i+1)=-100000000000;
%         elseif j>4*um*m+bm*p
%             bk(j,i+1)=-100000000000;
        end
    end
    yyy=eig(2.*H);
    Ik=[eye(um) zeros(um,um*(m-1))];
    xxx=quadprog(2.*H,-1.*G(:,i+1),-1.*Cu,-1.*bk(:,i+1));
    du(:,i)=Ik*xxx;

    if i==1
        u(:,i)=u(:,i)+du(:,i);%u(-1)=0
    else
        u(:,i)=u(:,i-1)+du(:,i);
    end   
    x(:,i+1)=Ad*x(:,i)+Bu*u(:,i)+Bd*road4(:,i);
    %***********被动比较***********%
    xp(:,i+1)=Ad*xp(:,i)+Bd*road4(:,i);%+Bu*u(:,i);
    yp(:,i)=Cc*xp(:,i);
end

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
plot(tout,x(19,:));
subplot(2,2,2)
plot(tout,x(20,:));
subplot(2,2,3)
plot(tout,x(21,:));
subplot(2,2,4)
plot(tout,x(22,:));