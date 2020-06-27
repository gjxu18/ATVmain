function [P,P12,P11]=transfer_uncontrollable_matrix_continuous_time(A,B,C,D,k,R,Q)

[SX,SY]=size(A);
state=SX;statey=SY;
[bx,by]=size(B);
rdd=sum(k);

A11=A(1:rdd,1:rdd);  
 Auc=A(rdd+1:state,rdd+1:state);
 A12=A(1:rdd,rdd+1:state);A21=A(rdd+1:state,1:rdd);A22=A(rdd+1:state,rdd+1:state);
B1=B(1:rdd,1:by);B2=B(rdd+1,1:by);

Rc=R(1:by,1:by);
C1=C(1:14,1:rdd);C2=C(1:14,rdd+1:state);
D1=D(1:rdd,1:4);  %D2=D(rdd+1:state,1:4);
%% 没有C
% Q11=Q(1:rdd,1:rdd); 
% Q22=Q(rdd+1:state,rdd+1:state);
% Q12=Q(1:rdd,rdd+1:state);
% Q21=Q(rdd+1:state,1:rdd);
% [P11,L,G]=care(A11,B1,Q11,Rc);
% % [jk,P11,rt]=lqr(A11,B1,Q11,R);
% P12=sylvester(A11'-P11*B1*Rc^(-1)*B1',A22,-Q12-P11*D1);
% P22=sylvester(A22',A22,-P12'*D1-D1'*P12+P12'*B1*R^(-1)*B1'*P12-Q22);
%% 论文的 有C 
[P11,L,G]=care(A11,B1,C1'*Q*C1,Rc);
P12=sylvester(A11'-P11*B1*Rc^(-1)*B1',A22,-C1'*Q*C2-P11*D1);
P22=sylvester(A22',A22,-P12'*D1-D1'*P12+P12'*B1*R^(-1)*B1'*P12-C2'*Q*C2);
%%
P21=P12';
%% calculate “P”
% Pnc=zeros(state-rdd);
% P=[P11,P12;
%     P21,Pnc];
P=[P11,P12;
    P21,P22];
