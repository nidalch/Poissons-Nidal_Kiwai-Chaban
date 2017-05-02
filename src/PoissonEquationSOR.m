% Poisson's Equation on a rectangle
% Nidal Kiwai Chaban
% SOR Method
clc
clear all
close all
n=input('Input number of nodes n for an n x n mesh: ');
x=n;
y=n;
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
u(1,:)=(((by-ay).^2).*cos(pi.*ay/by))+((xd(:)-ax)/(bx-ax)).*(((ay.*((by-ay).^2))-(((by-ay).^2.*cos(pi*ay/by)))));
u(y,:)=by;
        
tic        
while max(max(err(:)))>=1e-8
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
err=abs((uold-unew))./(unew);
end
el=toc;
%Inital BC


[X,Y]=meshgrid(xd,yd);
figure(1)
surf(X,Y,u) %3D Plot
xlabel('X domain')
ylabel('Y domain')
title('Successive Over Relaxation Solving of Poissons equation - Nidal Kiwai Chaban')
figure(2)
contourf(u) %2D Plot
xlabel('X domain (Nodes)')
ylabel('Y domain (Nodes)')
title('Successive Over Relaxation Solving of Poissons equation - Nidal Kiwai Chaban')

erb=max(max(err));
%Reporting

fprintf('Succesive-Over-Relaxation Method %d')
fprintf('\nMesh size: %d',x)
fprintf(' x %d',y)
fprintf('\nIterations: %d',k)
fprintf('\nElapsed time: %f',el);
fprintf('\nBiggest error value between u(i) and u(i-1): %.4ef',erb);