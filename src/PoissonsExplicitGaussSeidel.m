% Poisson's Equation on a rectangle
% Nidal Kiwai Chaban
% Gauss-Seidel Method
clc
clear all
close all
x=100;
y=100;
u=zeros(x,y);
bx=pi;
ax=-pi;
by=pi;
ay=-pi;
xd=linspace(ax,bx,x);
yd=linspace(ay,by,y);
ni=100
h=ax/x;
for k=1:ni
for i=2:x-1
    for j=2:y-1
        %F(i,j)=0;
        F(i,j)=cos((pi/2)*(2*((xd(i)-ax)/(bx-ax))+1))*sin(pi*(yd(j)-ay)/(by-ay));
        u(i,j)=(1/4)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)+(h.^2)*F(i,j));
        
        u(1,j)=((by-yd(j)).^2)*cos(pi*yd(j)/by);
        u(x,j)=yd(j)*(by-yd(j)).^2;
        u(i,y)=yd(j);
        u(i,1)=((by-ay).^2)*cos(pi*ay/by)+((xd(i)-ax)/(bx-ax))*((ay*((by-ay).^2)-((by-ay).^2)*cos(pi*ay/by)));
       
    end
end
end
%Inital BC
 u(1,1)=((by-yd(1)).^2)*cos(pi*yd(1)/by)
 u(x,1)=yd(1)*(by-yd(1)).^2;
 u(x,1)=yd(1)*(by-yd(1)).^2;


[X,Y]=meshgrid(xd,yd);
figure(1)
surf(X,Y,u) %3D Plot
xlabel('X domain')
ylabel('Y domain')
zlabel('Position')
title('Gauss-Seidel Solving of Poissons eqn')
figure(2)
contourf(u) %2D Plot
fprintf('Iterations: %f',ni)