function [Kmpc,Sx,I,Sd]=Model_Predictive_Control(Ad,Bu,Bd,Cc,Ts,p,m,rho_y,rho_u)
% p=20;m=10;rho_y=ones(1,3);rho_u=ones(1,4);
%*********增量状态空间模型********%
[xm,xm]=size(Ad);[xm,um]=size(Bu);[xm,dm]=size(Bd);
[ym,xm]=size(Cc);
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
Ik=[eye(um) zeros(um,um*(m-1))];
Kmpc1=pinv(Su'*(Rho_y'*Rho_y)*Su+(Rho_u'*Rho_u));
Kmpc2=(Su')*(Rho_y')*Rho_y;
Kmpc=Ik*Kmpc1*Kmpc2;








