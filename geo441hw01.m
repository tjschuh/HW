function geo441hw01(n)
% GEO441(n)
%
% code for HW 1
%
% n       what problem do you want to run
%         1 for problem 1a
%         2 for problem 1b
%         3 for problem 2
%
% Originally written by tschuh-at-princeton.edu, 02/02/2022
% Last modified by tschuh-at-princeton.edu, 02/07/2022

switch n
    case 1 % homogeneous, Dirichlet BCs
      % wave speed
      c = 1;

      % grid size and timestep
      dx = 0.1; dt = dx/c;

      % string length and max time
      xmax = 100; tmax = 2500;

      % create actual string
      x = [0:dx:xmax];
      
      % only ever need 3 rows
      % previous timestep --> current timestep --> future timestep
      u = zeros(3,xmax/dx+1);

      % define each row of displacement grid
      old = 1; cur = 2; new = 3;
      
      % IC: compute displacement values for t = 0
      for i=1:size(u,2)
          u(cur,i) = exp(-0.1*(((i-1)/10) - 50)^2);
      end

      % displacements for t = -1 are equal to t = 0 values
      % need this for 1st timestep
      u(old,:) = u(cur,:);
      
      % Dirichlet BCs (fixed ends)
      u(:,1) = 0; u(:,end) = 0;

      % plot t = 0 displacements
      f=figure;
      f.Visible = 'off';
      counter = 1;
      pint = 25;
      frate = 6;      
      plot(x,u(cur,:),'LineWidth',2)
      ylim([-1 1])
      M(counter) = getframe;
      counter = counter + 1;
      
      % compute future times aka "new" row by using equation from class
      for j=1:tmax
          for k=2:size(u,2)-1
              u(new,k) = ((c*dt/dx)^2)*(u(cur,k+1) - 2*u(cur,k) + u(cur,k-1)) + 2*u(cur,k) - u(old,k);
          end
          % after computing new displacements
          % current disps become old disps and new disps become cur disps
          u(old,:) = u(cur,:); u(cur,:) = u(new,:);
          % make movie and only plot every pint frame
          if mod(j,pint+1) == 1
              plot(x,u(cur,:),'LineWidth',2)
              ylim([-1 1])
              M(counter) = getframe;
              counter = counter + 1;
          end
      end

      % play movie
      f.Visible = 'on';
      movie(M,1,frate);

    case 2 % homogeneous, Neumann BCs, velocity-stress formulation
      % kappa, rho, and wave speed
      k = 1; p = 1; c = 1;

      % grid size and timestep
      dx = 0.1; dt = dx/c;

      % string length and max time
      xmax = 100; tmax = 100;

      % create actual string
      x = [0:dx:xmax];

      % allocate displacment array for first timestep only
      u = zeros(1,xmax/dx+1);

      % IC
      % compute displacement values for first row aka t = 0
      for i=1:size(u,2)
          u(1,i) = exp(-0.1*(((i-1)/10) - 50)^2);
      end

      % allocate velocity grid
      % we know that v(t=0) = 0
      v = zeros(tmax/dt,xmax/dx+1);

      % allocate stress array
      T = zeros(tmax/dt,xmax/dx+1);

      % Neumann BCs (stress-free ends)
      T(:,1) = 0;
      T(:,end) = 0;
      
      % now compute stress values for t = 0 using displacements at t=0
      for i=2:size(T,2)-1
          T(1,i) = (k/(2*dx))*(u(1,i+1) - u(1,i-1));
      end

      keyboard
      
      % now iterate to compute v and T for all timesteps
      % this is correct, but what are BCs on v (v(:,1) = v(:,end) = ???)
      for i=1:size(v,1)-1
          for j=2:size(T,2)-1
              if i == 1
                  v(i+1,j) = (dt/(p*dx))*(T(i,j+1) - T(i,j-1)) + v(i,j);
                  T(i+1,j) = (k/(2*dx))*(v(i,j+1) - v(i,j-1)) + T(i,j);
              else
                  v(i+1,j) = (dt/(p*dx))*(T(i,j+1) - T(i,j-1)) + v(i-1,j);
                  T(i+1,j) = (k/(2*dx))*(v(i,j+1) - v(i,j-1)) + T(i-1,j);
              end
          end
      end

      %keyboard
      
      % turn plots into movie
      % this is slow, can play with frame rate and plotting interval
      f=figure;
      f.Visible = 'off';
      counter = 1;
      int = 10;
      frate = 1;
      for m = 1:int:size(v,1)
          plot(x,T(m,:),'LineWidth',2)
          ylim([-1 1])
          drawnow
          M(counter) = getframe;
          counter = counter + 1;
      end
      f.Visible = 'on';
      movie(M,1,frate);
                  
    case 3 % heterogeneous, Dirichlet BCs NEED TO CHANGE TO Dirichlet at x=0 and Neumann at x=100
      % wave speed for [0 60] and (60 100]
      c1 = 1; c2 = 2;

      % grid size and timestep
      % use c2 when computing dt b/c we want smallest possible dt
      dx = 0.1; dt = dx/c2;

      % string length and max time
      xmax = 100; tmax = 100;

      % create actual string
      x = [0:dx:xmax];
      
      % allocate displacment grid
      u = zeros(tmax/dt,xmax/dx+1);

      % IC
      % compute displacement values for first row aka t = 0
      for i=1:size(u,2)
          u(1,i) = exp(-0.1*(((i-1)/10) - 50)^2);
      end

      % Dirichlet BCs (fixed ends)
      u(:,1) = 0;
      u(:,end) = 0;

      % compute future times aka subsequent rows by using equation from class
      for j=1:size(u,1)-1 % go until the 2nd to last time row bc we need j+1 terms
          for k=2:size(u,2)-1 % only go from 2nd element to 2nd to last element bc ends are fixed
              % for row 2 (t = 1), j - 1 term is the same as the j term
              if j == 1
                  if k <= 600
                      u(j+1,k) = ((c1*dt/dx)^2)*(u(j,k+1) - 2*u(j,k) + u(j,k-1)) + 2*u(j,k) - u(j,k);
                  elseif k == 601
                      u(j+1,k) = (dt/dx)^2*[(c2^2)*u(j,k+1) + (c1^2)*(-2*u(j,k) + u(j,k-1))] + 2*u(j,k) - u(j,k);
                  elseif k == 602
                      u(j+1,k) = (dt/dx)^2*[(c2^2)*(u(j,k+1) - 2*u(j,k)) + (c1^2)*u(j,k-1)] + 2*u(j,k) - u(j,k);
                  else %k > 603
                      u(j+1,k) = ((c2*dt/dx)^2)*(u(j,k+1) - 2*u(j,k) + u(j,k-1)) + 2*u(j,k) - u(j,k);
                  end
              % for every following timestep, actually use j - 1 term    
              else
                  if k <= 600
                      u(j+1,k) = ((c1*dt/dx)^2)*(u(j,k+1) - 2*u(j,k) + u(j,k-1)) + 2*u(j,k) - u(j-1,k);
                  elseif k == 601
                      u(j+1,k) = (dt/dx)^2*[(c2^2)*u(j,k+1) + (c1^2)*(-2*u(j,k) + u(j,k-1))] + 2*u(j,k) - u(j-1,k);
                  elseif k == 602
                      u(j+1,k) = (dt/dx)^2*[(c2^2)*(u(j,k+1) - 2*u(j,k)) + (c1^2)*u(j,k-1)] + 2*u(j,k) - u(j-1,k);
                  else %k > 603
                      u(j+1,k) = ((c2*dt/dx)^2)*(u(j,k+1) - 2*u(j,k) + u(j,k-1)) + 2*u(j,k) - u(j-1,k);
                  end
              end
          end
      end

      % turn plots into movie
      % this is slow, can play with frame rate and plotting interval
      f=figure;
      f.Visible = 'off';
      counter = 1;
      int = 10;
      frate = 10;
      for m = 1:int:size(u,1)
          plot(x,u(m,:),'LineWidth',2)
          xline(60,'--','LineWidth',2)
          text(27,0.8,'c = 1')
          text(77,0.8,'c = 2')
          ylim([-1 1])
          drawnow
          M(counter) = getframe;
          counter = counter + 1;
      end
      f.Visible = 'on';
      movie(M,1,frate);
end