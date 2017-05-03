% Poisson's Equation on a rectangle
% Nidal Kiwai Chaban
% Gauss-Seidel Method
clc
clear all
n=input('Input number of nodes n for an n x n mesh: ');
p=input('Input number of iterations: ');
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

%Boundary conditions
 u(:,1)=((by-yd(:)).^2).*cos(pi.*yd(:)/by);
 u(:,x)=yd(:).*(by-yd(:)).^2;
 u(1,:)=(((by-ay).^2).*cos(pi.*ay/by))+((xd(:)-ax)/(bx-ax)).*((ay.*((by-ay).^2)-((by-ay).^2).*cos(pi*ay/by)));
 u(y,:)=by;

tic 
%while max(max(err(:)))>=1e-6     %Run for optimal iterations
for q=1:p                        %Run for fixed iterations
    k=k+1;
    uold=u;
    for i=2:x-1
    for j=2:y-1
        %F(i,j)=0;
        F(i,j)=cos((pi/2)*(2*((xd(i)-ax)/(bx-ax))+1))*sin(pi*(yd(j)-ay)/(by-ay));
        u(i,j)=(1/4)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)+(h.^2)*F(i,j));
        
    end

    end
unew=u;
err=abs((uold-unew)./unew);
end
el=toc;


[X,Y]=meshgrid(xd,yd);

figure(1)             %Surface 3D Plot
surf(X,Y,u)
xlabel('X domain')
ylabel('Y domain')
zlabel('U Position') 
title(['Gauss-Seidel Solving of Poissons equation with '  num2str(n) ' x ' num2str(n)  ' mesh - Nidal Kiwai Chaban '])
colorbar

figure(2)
contourf(u) %2D Plot
xlabel('X domain (Nodes)')
ylabel('Y domain (Nodes)')
title(['Gauss-Seidel Solving of Poissons equation with '  num2str(n) ' x ' num2str(n)  ' mesh - Nidal Kiwai Chaban '])

erb=max(max(err));

%Reporting
%
fprintf('\nGauss-Seidel Method for Poissons Equation %d')
fprintf('\nMesh size: %d',x)
fprintf(' x %d',y)
fprintf('\nOptimal Iterations: %d',k)
fprintf('\nFixed Iterations: %d',p)
fprintf('\nElapsed time: %f',el);
fprintf('\nBiggest error value between u(i) and u(i-1): %10.4e',erb);