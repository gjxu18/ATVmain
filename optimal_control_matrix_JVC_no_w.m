function [Q,R]=optimal_control_matrix_JVC(parameter,A,B,SX,SY,ro)
%% 由cost function 计算Q R or N
syms x1 x2 x3 x4 x5 x6 x7 x8 x9...
    x10 x11 x12 x13 x14 x15 x16 x17 x18...
    x19 x20 x21 x22 x23 x24 x25 x26 x27...
    x28 x29 x30 x31 x32 x33 x34 x35 x36...
     w1 w2 w3 w4 u1 u2 u3 u4...
     a b c d Bf Br...
     ro1 ro2 ro3 ro4 ro5 ro6 ro7...%1垂向加速度、2tire deflections、3suspension deflections、4角度积分项、5u、6 7角度项;
     ro8 ro9 ro10...                     %8 9角加速度项 10垂向位置
     dx1 dx2 dx3 dx4
% ro1=ro(1); ro2=ro(2); ro3=ro(3); ro4=ro(4); ro5=ro(5);ro6=ro(6);ro7=ro(7);ro8=ro(8);ro9=ro(9);
a=parameter(5); b=parameter(6); Bf=parameter(7); Br=parameter(8);
zb1=x4-a*x5+1/2*Bf*x6;
zb2=x4-a*x5-1/2*Br*x6;
zb3=x4+b*x5+1/2*Br*x6;
zb4=x4+b*x5-1/2*Br*x6;
% dx1=(25*x15)/2 - (25*x8)/2 - (25*x9)/2 - (25*x10)/2 - 225*x11 - 225*x12 - 225*x13 - 225*x14 - (25*x7)/2 + (25*x16)/2 + (25*x17)/2 + (25*x18)/2 + 225*x19 + 225*x20 + 225*x21 + 225*x22;
% dx2=(2200*x7)/73 + (2200*x8)/73 - (2200*x9)/73 - (2200*x10)/73 + (39600*x11)/73 + (39600*x12)/73 - (39600*x13)/73 - (39600*x14)/73 - (2200*x15)/73 - (2200*x16)/73 + (2200*x17)/73 + (2200*x18)/73 - (39600*x19)/73 - (39600*x20)/73 + (39600*x21)/73 + (39600*x22)/73;
% dx3=- 40*x7 + 40*x8 - 40*x9 + 40*x10 - 720*x11 + 720*x12 - 720*x13 + 720*x14 + 40*x15 - 40*x16 + 40*x17 - 40*x18 + 720*x19 - 720*x20 + 720*x21 - 720*x22;

%% cost function 
% J=ro1*dx1^2+...
%     ro2*(w1-x19)^2+ro2*(w2-x20)^2+ro2*(w3-x21)^2+ro2*(w4-x22)^2+...
%     ro3*(x19-zb1)^2+ro3*(x20-zb2)^2+ro3*(x21-zb3)^2+ro3*(x22-zb4)^2+...
%     ro4*x31^2+ro4*x32^2+...
%     ro5*u1^2+ro5*u2^2+ro5*u3^2+ro5*u4^2+...
%     ro6*x5^2+ro7*x6^2+...
%     ro8*dx2^2+ro9*dx3^2+...
%     ro10*x4;
%% x不带w的
% J=ro(1)*dx1^2+ro(2)*dx2^2+ro(3)*dx3^2+ro(4)*x4^2+ro(5)*x5^2+ro(6)*x6^2+...
%     ro(7)*(w1-x19)^2+ro(8)*(w2-x20)^2+ro(9)*(w3-x21)^2+ro(10)*(w4-x22)^2+...
%     ro(11)*(x19-zb1)^2+ro(12)*(x20-zb2)^2+ro(13)*(x21-zb3)^2+ro(14)*(x22-zb4)^2+...
%     ro(15)*u1^2+ro(15)*u2^2+ro(15)*u3^2+ro(15)*u4^2;
%% x带w的
J=ro(1)*dx1^2+ro(2)*dx2^2+ro(3)*dx3^2+ro(4)*x4^2+ro(5)*x5^2+ro(6)*x6^2+...
    ro(7)*(x33-x19)^2+ro(8)*(x34-x20)^2+ro(9)*(x35-x21)^2+ro(10)*(x36-x22)^2+...
    ro(11)*(x19-zb1)^2+ro(12)*(x20-zb2)^2+ro(13)*(x21-zb3)^2+ro(14)*(x22-zb4)^2+...
    ro(15)*u1^2+ro(15)*u2^2+ro(15)*u3^2+ro(15)*u4^2;


% cf1=collect(J,[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32 x33 x34 x35 x36 u1 u2 u3 u4]);


 
%% 
xq=[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32 x33 x34 x35 x36];
uw=[u1 u2 u3 u4 w1 w2 w3 w4;];
%% Q
Q=zeros(SX,SY);
for i=1:SX
    cfi=diff(J,xq(1,i));
    for j=1:SY
        if j>=i
            Q(i,j)=1/2*diff(cfi,xq(1,j))+A(1,i)*A(1,j)*ro(1)+A(2,i)*A(2,j)*ro(2)+A(3,i)*A(3,j)*ro(3); 
        end
    end
end
% for i=1:SX
%     for j=1:SY
%         if j>=i
%             Q(i,j)=A(1,i)*A(1,j)*ro1;
%         end
%     end
% end                                   %%定义ro1的Q
for i=1:SX
    for j=1:SY
        if j>i
            Q(j,i)=Q(i,j);
        end
    end
end
%% R
[bm,bn]=size(B);
R=zeros(bn,bn);
for i=1:bn
    cfu=diff(J,uw(1,i));
    for j=1:bn
        if j>=i
            R(i,j)=1/2*diff(cfu,uw(1,j));%+A(1,i)*A(1,j)*ro1+A(2,i)*A(2,j)*ro8+A(3,i)*A(3,j)*ro9;  %costfunction系数+定义rol的Q+角加速度项ro8 ro9
        end
    end
end
for i=1:bn
    for j=1:bn
        if j>i
            R(j,i)=R(i,j);
        end
    end
end
%%
% Q(4,4)=Q(4,4)+4*ro3;
% Q(4,5)=Q(4,5)+2*(b-a)*ro3;
% Q(4,6)=Q(4,6)+1/2*(Bf*ro3 - Br*ro3);
% Q(4,19)=Q(4,19)-ro3;
% Q(4,20)=Q(4,20)-ro3;
% Q(4,21)=Q(4,21)-ro3;
% Q(4,22)=Q(4,22)-ro3;
% 
% Q(5,5)=Q(5,5)+(2*ro3*a^2+2*ro3*b^2);
% Q(5,19)=Q(5,19)+a*ro3;
% Q(5,6)=Q(5,6)+1/2*(Br*a*ro3 - Bf*a*ro3);
% Q(5,20)=Q(5,20)+a*ro3;
% Q(5,21)=Q(5,21)-b*ro3;
% Q(5,22)=Q(5,22)-b*ro3;
% 
% Q(6,6)=Q(6,6)+((ro3*Bf^2)/4 + (3*ro3*Br^2)/4);
% Q(6,19)=Q(6,19)-1/2*Bf*ro3;
% Q(6,20)=Q(6,20)+1/2*Br*ro3;
% Q(6,21)=Q(6,21)-1/2*Br*ro3;
% Q(6,22)=Q(6,22)+1/2*Br*ro3;
% %% 
% Q(5,5)=Q(5,5)+ro6;
% Q(6,6)=Q(6,6)+ro6;
% %%
% Q(19,19)=Q(19,19)+(ro2+ro3);
% Q(19,33)=Q(19,33)-ro2;
% Q(20,20)=Q(20,20)+(ro2+ro3);
% Q(20,34)=Q(20,34)-ro2;
% Q(21,21)=Q(21,21)+(ro2+ro3);
% Q(21,35)=Q(21,35)-ro2;
% Q(22,22)=Q(22,22)+(ro2+ro3);
% Q(22,36)=Q(22,36)-ro2;
% 
% Q(31,31)=Q(31,31)+ro4;
% Q(32,32)=Q(32,32)+ro4;
% 
% Q(33,33)=Q(33,33)+ro2;
% Q(34,34)=Q(34,34)+ro2;
% Q(35,35)=Q(35,35)+ro2;
% Q(36,36)=Q(36,36)+ro2;
% 
% for i=1:SX
%     for j=1:SY
%         if j>i
%             Q(j,i)=Q(i,j);
%         end
%     end
% end
%%
% [bm,bn]=size(B);
% R=zeros(bn,bn);
% for i=1:bn
%     for j=1:bn
%         R(i,j)=B(1,i)*B(1,j)*ro1;
%     end
% end
% R(1,1)=R(1,1)+ro5;
% R(2,2)=R(2,2)+ro5;
% R(3,3)=R(3,3)+ro5;
% R(4,4)=R(4,4)+ro5;
