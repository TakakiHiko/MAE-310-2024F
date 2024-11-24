clear;
clc;
ux1=0;
u0=1;
u=@(x) x.^5;
calu=1/11; %u平方在0到1的积分
calux=16/9;%u对x的偏导的平方在0到1的积分值
max_point=16;
erro=zeros(max_point/2,2);%create an arrary to storage error

for i=2:2:16
pp=1;
n_el=i;
n_int=10;
uh= driver(pp,n_el,n_int,ux1,u0);
h=1/n_el/pp;%Determine the interval between each adjacent point
x=0:h:1-h;
ur=u(x);

for j=1:i
%     erro(i/2,1)=erro(i/2,1)+;
end
end





