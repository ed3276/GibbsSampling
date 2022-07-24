tic
N=100;
x=linspace(-5,10,N+1);y=linspace(-5,10,N+1);
[X,Y]=meshgrid(x,y);
pa.u1=[0;0];pa.sgm1=[1.5,0.2;0.2,1.4];pa.rou1=pa.sgm1(1,2)/(sqrt(pa.sgm1(1,1))*sqrt(pa.sgm1(2,2)));
pa.u2=[2;3];pa.sgm2=[1.1,0.4;0.4,1.3];pa.rou2=pa.sgm2(1,2)/(sqrt(pa.sgm2(1,1))*sqrt(pa.sgm2(2,2)));
Z=zeros(N,N);
for j=1:length(y)
    for i=1:length(x)
        Z(j,i)=fxyv(pa,[x(i);y(j)]);
    end
end
subplot(2,2,1)
surf(X,Y,Z)
view(15,30)
axis([-5,10,-5,10,0,0.1]);
xlabel('x');ylabel('y');
hold on
load randGibbs.dat
subplot(2,2,2)
histogram2(randGibbs(:,1)',randGibbs(:,2)','Normalization','pdf','FaceColor','flat');
view(15,30)
axis([-5,10,-5,10,0,0.1]);
xlabel('x');ylabel('y');

subplot(2,2,3)
surf(X,Y,Z)
view(0,90)
axis([-5,10,-5,10,0,0.1]);
xlabel('x');ylabel('y');
hold on
subplot(2,2,4)
histogram2(randGibbs(:,1)',randGibbs(:,2)','Normalization','pdf','FaceColor','flat');
view(0,90)
axis([-5,10,-5,10,0,0.1]);
xlabel('x');ylabel('y');

function y=fxyv(pa,x)
    u1=pa.u1;sgm1=pa.sgm1;
    u2=pa.u2;sgm2=pa.sgm2;
    y=1/2*fxy(u1,sgm1,x)+1/2*fxy(u2,sgm2,x);
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
end
function y=fxy(u,s,x)
    y=1/(2*pi*sqrt(det(s)))*exp(-(x-u)'/s*(x-u)/2);
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
end
function y=f(u,s,x)
    y=1/(sqrt(2*pi)*s)*exp(-(x-u)^2/(s^2*2));
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
end
function y=fxv(pa,x)
    u1=pa.u1;sgm1=pa.sgm1;
    u2=pa.u2;sgm2=pa.sgm2;
    y=1/2*f(u1(1),sqrt(sgm1(1,1)),x)+1/2*f(u2(1),sqrt(sgm2(1,1)),x);
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
end
function y=fyv(pa,x)
    u1=pa.u1;sgm1=pa.sgm1;
    u2=pa.u2;sgm2=pa.sgm2;
    y=1/2*f(u1(2),sqrt(sgm1(2,2)),x)+1/2*f(u2(2),sqrt(sgm2(2,2)),x);
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
end