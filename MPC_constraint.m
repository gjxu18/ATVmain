function [Kg,Sx,I,Sd,Sxb,Sub,Sdb,Ib,H,Cu]=MPC_constraint(Ad,Bu,Bd,Cc,Cb,Ts,p,m,rho_y,rho_u)
% p=20;m=10;rho_y=ones(1,3);rho_u=ones(1,4);
%*********增量状态空间模型********%
[xm,xm]=size(Ad);[xm,um]=size(Bu);[xm,dm]=size(Bd);
[ym,xm]=size(Cc);[bm,xm]=size(Cb);
%*********预测方程*********%
Sx=zeros(ym*p,xm);
I=zeros(ym*p,ym);                   
Sd=zeros(ym*p,dm);
Su=zeros(ym*p,um*m);
for i=1:p
    Sxj=zeros(ym,xm);
    Sdj=zeros(ym,dm);
    for j=1:i
        Sxj=Sxj+Cc*Ad^j;
        Sdj=Sdj+Cc*Ad^(j-1)*Bd;
    end
    Sx((i-1)*ym+1:ym*i,:)=Sxj;
    I((i-1)*ym+1:ym*i,:)=eye(ym);
    Sd((i-1)*ym+1:ym*i,:)=Sdj;
    for k=1:m
        if i>=k
            Suj=zeros(ym,um);
            for j=1:i-k+1
                Suj=Suj+Cc*Ad^(j-1)*Bu;
            end
            Su((i-1)*ym+1:ym*i,(k-1)*um+1:um*k)=Suj;
        end
    end
end
%***********控制约束**************
T=eye(m*um);
L=zeros(m*um,m*um);
for i=1:m
    for h=1:m
        if i>=h
            IL=eye(um);
            L((i-1)*um+1:um*i,(h-1)*um+1:um*h)=IL;
        end
    end
end
%****************输出约束**************
Sxb=zeros(bm*p,xm);
Ib=zeros(bm*p,bm);                   
Sdb=zeros(bm*p,dm);
Sub=zeros(bm*p,um*m);
for i=1:p
    Sxbj=zeros(bm,xm);
    Sdbj=zeros(bm,dm);
    for j=1:i
        Sxbj=Sxbj+Cb*Ad^j;
        Sdbj=Sdbj+Cb*Ad^(j-1)*Bd;
    end
    Sxb((i-1)*ym+1:ym*i,:)=Sxbj;
    Ib((i-1)*ym+1:ym*i,:)=eye(bm);
    Sdb((i-1)*ym+1:ym*i,:)=Sdbj;
    
    for k=1:m
        if i>=k
            Subj=zeros(ym,um);
            for j=1:i-k+1
                Subj=Subj+Cb*Ad^(j-1)*Bu;
            end
            Sub((i-1)*ym+1:ym*i,(k-1)*um+1:um*k)=Subj;
        end
    end
end
              
%*********加权因子&参考输出*******%
rho_y=rho_y;
rho_u=rho_u;
Rho_y=zeros(1,ym*p);
Rho_u=zeros(1,um*m);
for i=1:p
    Rho_y((i-1)*ym+1:ym*i)=rho_y;     %可改
end
for i=1:m
    Rho_u((i-1)*um+1:um*i)=rho_u;     %可改
end
%***********开环优化************%
H=Su'*(Rho_y'*Rho_y)*Su+(Rho_u'*Rho_u);
Cu=[-T' T' -L' L' -Sub' Sub']';
Cu=[-T' T' -L' L']';
Kg=2.*Su'*(Rho_y'*Rho_y);
% H=(H+H')/2;

Ik=[eye(um) zeros(um,um*(m-1))];
Kmpc1=pinv(((Su')*(Rho_y')*Rho_y*Su+(Rho_u')*Rho_u),1e-8);
Kmpc2=(Su')*(Rho_y')*Rho_y;
Kmpc=Ik*Kmpc1*Kmpc2;








