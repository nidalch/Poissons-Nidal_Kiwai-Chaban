% Poisson's Equation on a rectangle
% Nidal Kiwai Chaban
% Gauss-Seidel Method
clc
clear all
close all
x=50;
y=50;
u=zeros(x,y);
bx=pi;
ax=-pi;
by=pi;
ay=-pi;
xd=linspace(ax,bx,x);
yd=linspace(ay,by,y);
ni=100
h=ax/x;

 u(1,:)=((by-yd(:)).^2).*cos(pi.*yd(:)/by);
 u(x,:)=yd(:).*(by-yd(:)).^2;
 u(:,y)=yd(:);
 u(:,1)=(((by-ay).^2).*cos(pi.*ay/by))+((xd(:)-ax)/(bx-ax)).*((ay.*((by-ay).^2)-((by-ay).^2).*cos(pi*ay/by)));

        
    for k=1:ni
    for i=2:x-1
    for j=2:y-1
        %F(i,j)=0;
        F(i,j)=cos((pi/2)*(2*((xd(i)-ax)/(bx-ax))+1))*sin(pi*(yd(j)-ay)/(by-ay));
        u(i,j)=(1/4)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)+(h.^2)*F(i,j));
        
    end
    u(i,y)=yd(i);
end
end
%Inital BC


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