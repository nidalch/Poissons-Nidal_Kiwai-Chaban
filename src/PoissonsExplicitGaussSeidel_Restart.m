%% Restarting
matfile = 'PoissonEquationSolution';   % should match that in test_checkpoint.m
load(matfile);        % retrieve data from matfile
k1 = k+1;       % iter is the last time test_checkpoint issued
                      % a save; we start computing on the next step

%% Iterative Loop 
tic; 
while max(err(:))>=1e-6     %Run for optimal iterations
%for q=1:p                        %Run for fixed iterations
    k=k+1;
    
        
    % If you want to test the restart script, use the function pause(1) to slow down the while loop.
    % This will slow down the while loop to 1 sec per iteration so that ctrl + C can used be to
    % "kill" the code to simulate a computer crash. From there, use the restart script to restart the loop.  
   %  pause(.01)
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
%             F(i,j)=cos((pi/2)*(2*((xd(i)-ax)/(bx-ax))+1))*sin(pi*(yd(j)-ay)/(by-ay));
            u(i,j)=(1/4)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)+(h.^2)*F(i,j));
            u(x,j)=(1/4)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j+1)+(h.^2)*F(i,j)); %Neumann BC
        end

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
title(['Gauss-Seidel Solving of Poissons equation with '  num2str(n) ' x ' num2str(n)  ' mesh - Nidal Kiwai Chaban '],'FontSize',9)

figure(2)
contourf(xd,yd,u) %2D Plot
xlabel('X domain')
ylabel('Y domain')
title(['Gauss-Seidel Solving of Poissons equation with '  num2str(n) ' x ' num2str(n)  ' mesh - Nidal Kiwai Chaban '],'FontSize',9)

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