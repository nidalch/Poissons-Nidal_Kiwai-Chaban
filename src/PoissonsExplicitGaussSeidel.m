% Poisson's Equation on a rectangle
% Nidal Kiwai Chaban
% Gauss-Seidel Method
clc
clear all
close all
x=40;
y=40;
u=zeros(x,y);
bx=pi;
ax=-pi;
by=pi;
ay=-pi;
ni=2000
h=ax/ni;
for k=1:ni
for i=2:x-1
    for j=2:y-1
%         F(i,j)=
        u(i,j)=(1/4)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1));
        u(:,x)=bx;
        u(:,1)=ax;
        u(1,:)=by;
        u(y,:)=ay;
    end
end
end
xp=1:x;
yp=1:y;
[X,Y]=meshgrid(xp,yp);
subplot(1,2,1)
surf(X,Y,u) %3D Plot
subplot(1,2,2)
contourf(u) %2D Plot
fprintf('Iterations: %f',ni)