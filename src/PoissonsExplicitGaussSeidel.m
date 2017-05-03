%% Poisson's Equation on a rectangle
% Nidal Kiwai Chaban
% Gauss-Seidel Method
clc
clear all

%% Define Initial Values
n=input('Input number of nodes n for an n x n mesh: ');
%p=input('Input number of iterations: ');
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
freq=10;
err=1;
h=abs(ax)/x;

%% Boundary conditions
 u(:,1)=((by-yd(:)).^2).*cos(pi.*yd(:)/by);
 u(:,x)=yd(:).*(by-yd(:)).^2;
 u(1,:)=(((by-ay).^2).*cos(pi.*ay/by))+((xd(:)-ax)/(bx-ax)).*((ay.*((by-ay).^2)-((by-ay).^2).*cos(pi*ay/by)));
 %u(y,:)=by;
 F=cos((pi/2)*(2*((xd-ax)/(bx-ax))+1))'*sin(pi.*(yd-ay)/(by-ay));
%F=zeros(x,y);
%% Checkpointing
% Sometimes files take a long time to run to completion. As a result, sometimes they crash due to a variety of reasons: power failure, walltime limit, scheduled shutdown, etc. 
% Checkpoint/Restarting has long been a common technique to tackle this issue. Checkpointing/Restarting essentially means saving data to disk periodically so that, if need be, 
% you can restart the job from the point at which your data was last saved. 

% Before the start of the iteration loop, "check-in" each variable
% that should be checkpointed in the event of restarting the job

matfile = 'PoissonEquationSolution.mat';     % mandatory; name of checkpoint mat-file
s = struct();                                % mandatory; create struct for checkpointing
s = chkin(s,{'k'});                       % mandatory; iter is iteration loop index
s = chkin(s,{'freq'});                  % mandatory; frequency is checkpointing period 
                                             % i.e., how often to perform a save

% continue until all variables are checked in. Note that you are only
% checking in the variables, they don't need to have been already defined

chkNames = fieldnames(s);    % the full list of variables to checkpoint
nNames = length(chkNames);   % number of variables in list

%% Iterative Loop 
tic; 
while max(err(:))>=1e-6     %Run for optimal iterations
%for q=1:p                        %Run for fixed iterations
    k=k+1;
    
        
    % If you want to test the restart script, use the function pause(1) to slow down the while loop.
    % This will slow down the while loop to 1 sec per iteration so that ctrl + C can used be to
    % "kill" the code to simulate a computer crash. From there, use the restart script to restart the loop.  
    % pause(.001)
    if mod(k, freq) == 0 % If statement, checkpoints periodically (determined by the frequency)
        chkpt                    % chkpt script performs checkpointing (save) every *frequency* iterations
        fprintf(1, ['Checkpointing frequency is every %2d iterations.' ...
          'Data updated at iteration %3d\n'], ...
          freq, k);      % Confirm after each checkpointing event 
    end
    
    uold=u;
    for i=2:x-1
        for j=2:y-1
            %F(i,j)=0;
           % F(i,j)=cos((pi/2)*(2*((xd(i)-ax)/(bx-ax))+1))*sin(pi*(yd(j)-ay)/(by-ay));
            u(i,j)=(1/4)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)+(h.^2)*F(i,j));
            u(x,j)=(1/4)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j+1)+(h.^2)*F(i,j)); %Neumann BC
        end
  %          u(x,j)=(1/4)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j+1)+(h.^2)*F(i,j)); %Neumann BC
    end
unew=u;
err=abs((uold-unew))./unew;
    fprintf(1, 'Completed iteration %d\n', k);
end
el=toc;

%% Plotting
[X,Y]=meshgrid(xd,yd);

figure(1)             %Surface 3D Plot
surf(X,Y,u,'EdgeColor','none')
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
l=mean(mean(u));

%% Reporting

fprintf('\nGauss-Seidel Method for Poissons Equation %d')
fprintf('\nMesh size: %d',x)
fprintf(' x %d',y)
fprintf('\nOptimal Iterations: %d',k)
%fprintf('\nFixed Iterations: %d',p)
fprintf('\nElapsed time: %f',el);
fprintf('\nBiggest error value between u(i) and u(i-1): %10.4e',erb);
fprintf('\nMean of u: %10.4e',l);