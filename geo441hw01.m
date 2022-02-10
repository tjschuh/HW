function geo441hw01(n,m)
% GEO441HW01(n,m)
%
% code for HW 1, makes and saves a movie of an oscillating string
% with varying boundary conditions and material properties
%
% INPUT:
%    
% n       what problem you want to run
%         1 --> problem 1
%         2 --> problem 2
% m       if problem 1, what part of a problem you want to run
%         1 --> part a
%         2 --> part b
%
% OUTPUT:
%
% .mp4 file of the movie
%    
% EXAMPLES:
%
% geo441hw01(1,1)
% geo441hw01(1,2)
% geo441hw01(2,[])   
%    
% Originally written by tschuh-at-princeton.edu, 02/02/2022
% Last modified by tschuh-at-princeton.edu, 02/09/2022

switch n
    case 1 % homogeneous case
      % material properties
      % wave speed c = sqrt(k/p)
      k = 1; p = 1;
      c = sqrt(k/p);

      % grid size and timestep
      dx = 0.1; dt = dx/c;

      % string length and max time
      xmax = 100; tmax = 2000;

      % create actual string
      x = [0:dx:xmax];
      
      % only ever need 3 rows
      % previous timestep --> current timestep --> future timestep
      u = zeros(3,xmax/dx+1);

      % allocate velocity grid
      v = zeros(3,xmax/dx+1);

      % allocate stress array
      T = zeros(3,xmax/dx+1);

      % define each row of displacement grid
      old = 1; cur = 2; new = 3;

      % IC: compute displacement values for t = 0
      for i=1:size(u,2)
          u(cur,i) = exp(-0.1*(((i-1)/10) - 50)^2);
      end

      % displacements for t = -1 are equal to t = 0 values
      % need this for 1st timestep
      u(old,:) = u(cur,:);

      % now compute stress values for t=0 using displacements at t=0
      for i=2:size(T,2)-1
          T(cur,i) = (k/(2*dx))*(u(cur,i+1) - u(cur,i-1));
      end
      T(old,:) = T(cur,:);

      % v is always zero at t=0 so we dont need to do any ICs there

      % plot t = 0 displacements
      f=figure;
      f.Visible = 'off';
      counter = 1;
      pint = 25;
      frate = 6;      
      plot(x,u(cur,:),'k','LineWidth',2)
      hold on
      plot(x,v(cur,:),'b','LineWidth',2)
      plot(x,T(cur,:),'r','LineWidth',2)
      ylim([-1 1])
      legend('Displacement','Velocity','Displacement')
      grid on
      M(counter) = getframe;
      hold off
      counter = counter + 1;

      switch m
          case 1 % Dirichlet BCs                  
            % Dirichlet BCs (fixed ends)
            u(:,1) = 0; u(:,end) = 0;
            v(:,1) = 0; v(:,end) = 0;
          
            % compute future times aka "new" row by using discretized equations from class
            for i=1:tmax
                for j=2:size(u,2)-1
                    u(new,j) = ((c*dt/dx)^2)*(u(cur,j+1) - 2*u(cur,j) + u(cur,j-1)) + 2*u(cur,j) - u(old,j);
                    v(new,j) = (dt/(p*dx))*(T(cur,j+1) - T(cur,j-1)) + v(old,j);
                    T(new,j) = ((k*dt)/dx)*(v(cur,j+1) - v(cur,j-1)) + T(old,j);
                end
                % u(ends) = v(ends) = 0, so T(ends) = T(adjacent)
                T(new,1) = T(new,2); T(new,end) = T(new,end-1);
                % after computing new displacements
                % current vals become old vals and new vals become cur vals
                u(old,:) = u(cur,:); u(cur,:) = u(new,:);
                v(old,:) = v(cur,:); v(cur,:) = v(new,:);
                T(old,:) = T(cur,:); T(cur,:) = T(new,:);
                % make movie and only plot every pint frame
                if mod(i,pint-1) == 0
                    plot(x,u(cur,:),'k','LineWidth',2)
                    hold on
                    plot(x,v(cur,:),'b','LineWidth',2)
                    plot(x,T(cur,:),'r','LineWidth',2)
                    ylim([-1 1])
                    title('Homogeneous, Dirichlet BCs')
                    legend('Displacement','Velocity','Stress')
                    grid on
                    M(counter) = getframe;
                    hold off
                    counter = counter + 1;
                end
            end

            % play and save movie
            f.Visible = 'on';
            v = VideoWriter('1a','MPEG-4');
            v.FrameRate = frate;
            open(v)
            movie(M,2,frate);
            writeVideo(v,M)
            close(v)

          case 2 % Neumann BCs
            % Neumann BCs (stress-free ends)
            T(:,1) = 0; T(:,end) = 0;
          
            % compute future times aka "new" row by using discretized equations from class
            for i=1:tmax
                for j=2:size(T,2)-1
                    u(new,j) = ((c*dt/dx)^2)*(u(cur,j+1) - 2*u(cur,j) + u(cur,j-1)) + 2*u(cur,j) - u(old,j);
                    v(new,j) = (dt/(p*dx))*(T(cur,j+1) - T(cur,j-1)) + v(old,j);
                    T(new,j) = ((k*dt)/dx)*(v(cur,j+1) - v(cur,j-1)) + T(old,j);
                end
                % BCs: T(ends) = 0, so u(ends) = u(adjacent) (same for v)
                u(new,1) = u(new,2); u(new,end) = u(new,end-1);
                v(new,1) = v(new,2); v(new,end) = v(new,end-1);
                % update timesteps
                u(old,:) = u(cur,:); u(cur,:) = u(new,:);
                v(old,:) = v(cur,:); v(cur,:) = v(new,:);
                T(old,:) = T(cur,:); T(cur,:) = T(new,:);
                if mod(i,pint-1) == 0
                    plot(x,u(cur,:),'k','LineWidth',2)
                    hold on
                    plot(x,v(cur,:),'b','LineWidth',2)
                    plot(x,T(cur,:),'r','LineWidth',2)
                    ylim([-1 1])
                    legend('Displacement','Velocity','Stress')
                    title('Homogeneous, Neumann BCs')                    
                    grid on
                    M(counter) = getframe;
                    hold off              
                    counter = counter + 1;
                end
            end

            % play movie
            f.Visible = 'on';
            v = VideoWriter('1b','MPEG-4');
            v.FrameRate = frate;
            open(v)
            movie(M,2,frate);
            writeVideo(v,M)
            close(v)
      end
    case 2 % heterogeneous case
      % boundary between 2 regions in hetereogeneous case
      split = 60;

      % grid size
      dx = 0.1;

      % string length and max time
      xmax = 100; tmax = 2000;

      % create actual string
      x = [0:dx:xmax];

      % define each row of displacement grid
      old = 1; cur = 2; new = 3;
      
      % allocate displacment grid
      u = zeros(3,xmax/dx+1);

      % allocate velocity grid
      v = zeros(3,xmax/dx+1);

      % allocate stress array
      T = zeros(3,xmax/dx+1);

      % Dirichlet BC at x = 0 and Neumann BC at x = 100
      u(:,1) = 0; T(:,end) = 0;

      % IC: compute displacement values for first row aka t = 0
      for i=1:size(u,2)
          u(cur,i) = exp(-0.1*(((i-1)/10) - 50)^2);
      end
      u(old,:) = u(cur,:);

      % material properties for [0 60] and (60 100]
      p = zeros(1,xmax/dx+1);
      p(1,1:(split/dx)) = 1;
      p(1,(split/dx+1):end) = 1;

      k = zeros(1,xmax/dx+1);
      k(1,1:(split/dx)) = 1;
      k(1,(split/dx+1):end) = 4;

      % wave speed c = sqrt(k/p)
      c = zeros(1,xmax/dx+1);
      c(1,1:(split/dx)) = sqrt(k(1:(split/dx))/p(1:(split/dx)));
      c(1,(split/dx+1):end) = sqrt(k((split/dx+1):end)/p((split/dx+1):end));      
      
      % we need smallest possible timestep so we need largest c value
      dt = dx/max(c);

      % now compute stress values for t=0 using displacements at t=0
      for i=2:size(T,2)-1
          T(cur,i) = (1/(2*dx))*(k(1,i+1)*u(cur,i+1) - k(1,i-1)*u(cur,i-1));
      end
      T(old,:) = T(cur,:);

      f=figure;
      f.Visible = 'off';
      counter = 1;
      pint = 25;
      frate = 6;      
      plot(x,u(cur,:),'k','LineWidth',2)
      hold on
      plot(x,v(cur,:),'b','LineWidth',2)
      plot(x,T(cur,:),'r','LineWidth',2)
      xline(split,'--','LineWidth',2)
      text(27,0.8,sprintf('c = %.f',c(1,split/dx)))
      text(77,0.8,sprintf('c = %.f',c(1,split/dx+1)))
      ylim([-1 1])
      legend({'Displacement','Velocity','Stress'},'Location','southeast')
      grid on
      M(counter) = getframe;
      hold off
      counter = counter + 1;

      % compute future times aka subsequent rows by using discretized equation from class with additional hetereogeneous term
      for i=1:tmax
          for j=2:size(u,2)-1 % only go from 2nd element to 2nd to last element because of BCs
              u(new,j) = ((dt/dx)^2)*[(c(1,j+1))^(2)*u(cur,j+1) - (c(1,j))^(2)*2*u(cur,j) + (c(1,j-1))^(2)*u(cur,j-1)] + 2*u(cur,j) - u(old,j);
              v(new,j) = (dt/dx)*((T(cur,j+1)/p(1,j+1)) - (T(cur,j-1)/p(1,j+1))) + v(old,j);
              T(new,j) = (dt/dx)*((k(1,j+1)*v(cur,j+1)) - (k(1,j-1)*v(cur,j-1))) + T(old,j);
          end
          % implement BCs
          % left-side of string u = v = 0, so T(x=0) = T(adjacent element)
          % right-side of string T = 0, so u(x=end) = u(adjacent element) (same for v)
          u(new,end) = u(new,end-1);
          v(new,end) = v(new,end-1);
          T(new,1) = T(new,2);
          % update timesteps
          u(old,:) = u(cur,:); u(cur,:) = u(new,:);
          v(old,:) = v(cur,:); v(cur,:) = v(new,:);
          T(old,:) = T(cur,:); T(cur,:) = T(new,:);
          if mod(i,pint-1) == 0
              plot(x,u(cur,:),'k','LineWidth',2)
              hold on
              plot(x,v(cur,:),'b','LineWidth',2)
              plot(x,T(cur,:),'r','LineWidth',2)
              xline(split,'--','LineWidth',2)
              text(27,0.8,sprintf('c = %.f',c(1,split/dx)))
              text(77,0.8,sprintf('c = %.f',c(1,split/dx+1)))
              ylim([-1 1])
                    title('Heterogeneous, Mixed BCs')              
              legend({'Displacement','Velocity','Stress'},'Location','southeast')
              grid on
              M(counter) = getframe;
              hold off              
              counter = counter + 1;
          end
      end

      % play movie
      f.Visible = 'on';
      v = VideoWriter('2','MPEG-4');
      v.FrameRate = frate;
      open(v)
      movie(M,2,frate);
      writeVideo(v,M)
      close(v)
end