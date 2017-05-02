% Poisson's Equation on a rectangle
% Nidal Kiwai Chaban
% SOR Method
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
k=0;
err=1;
h=ax/x;
B=1.5;

%Boundary conditions
 u(:,1)=((by-yd(:)).^2).*cos(pi.*yd(:)/by);
 u(:,x)=yd(:).*(by-yd(:)).^2;
 u(1,:)=(((by-ay).^2).*cos(pi.*ay/by))+((xd(:)-ax)/(bx-ax)).*((ay.*((by-ay).^2)-((by-ay).^2).*cos(pi*ay/by)));
 u(y,:)=by;
        
        
while max(err(:))>=1e-6
    k=k+1;
    uold=u;
    for i=2:x-1
    for j=2:y-1
        %F(i,j)=0;
        F(i,j)=cos((pi/2)*(2*((xd(i)-ax)/(bx-ax))+1))*sin(pi*(yd(j)-ay)/(by-ay));
        u(i,j)=(B.*(1/4)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)+(h.^2)*F(i,j)))+(1-B)*u(i,j);
        
    end
    end
unew=u;
err=abs((uold-unew)./unew);
end
%Inital BC


[X,Y]=meshgrid(xd,yd);
figure(1)
surf(X,Y,u) %3D Plot
xlabel('X domain')
ylabel('Y domain')
zlabel('Position')
title('SOR Solving of Poissons eqn')
figure(2)
contourf(u) %2D Plot
fprintf('Iterations: %f',k)