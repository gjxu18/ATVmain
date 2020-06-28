clc;
clear;
%% Êµ¼Ê
x1=(0:0.01:1.5);        %1.5>a+b T=6
x2=(1.51:0.01:4.5);     
x3=(4.51:0.01:6.0);
y1=tan(8/180*pi)*x1;    %8deg<(a+b)/0.2
y2=(3.0-x2)*tan(8/180*pi);
y3=-tan(8/180*pi)*(6.0-x3);
x=[x1';x2';x3'];
y=[y1';y2';y3'];
plot(x,y);
hold on
%% sin
a_max=tan(8/180*pi)/(2*pi/6); %0.1342
y_sin=a_max*sin(2*pi/6*x);  %y=0.1342sin(pi/3*x) T=6

plot(x,y_sin);


