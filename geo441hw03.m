function geo441hw03(n)
% GEO441HW03(n)
%    
% Given a problem number (1,2,3), use the respective finite difference
% scheme (forward, backward, Crank-Nicolson) to solve the heat equation
% and plot the temperature profile
%
% INPUT:
%
% n       what problem you want to run
%         1 --> problem 1
%         2 --> problem 2
%         3 --> problem 3
%
% OUTPUT:
%
% .mp4 file of the movie
%
% Originally written by tschuh-at-princeton.edu, 02/17/2022
% Last modified by tschuh-at-princeton.edu, 02/23/2022

% grid size and rod length
dx = 1; L = 100;

% create actual rod, N = (L/dx) + 1
x = [0:dx:L];

switch n
    case 1 % problem 1 (heterogeneous)
      fprintf('Working on Problem 1\n')

      % define temperature array
      T = zeros(2,length(x));
      old = 1; new = 2;

      % ICs
      T(old,L/2+1) = 1;

      % define thermal conductivity array
      k = zeros(1,length(x));
      k(:) = 1;
      % define density array
      p = zeros(1,length(x));
      p(:) = 1;
      % define specific heat array
      cp = zeros(1,length(x));
      cp(:) = 1;

      % define D s.t. it is as large as possible
      % so that dt is as small as possible
      D = max(k)/(min(p)*min(cp));

      % timestep and tmax
      % coeff > 0.5 crashes
      coeff = 0.4;
      dt = coeff*dx*dx/D;
      tmax = 150;

      % BCs
      T(:,1) = 0; T(:,length(x)) = 0;

      % initial plotting
      f=figure;
      f.Visible = 'off';
      counter = 1;
      pint = 5;
      plays = 1;
      frate = 4;
      plot(x,T(old,:),'k','LineWidth',2)
      xlabel('rod')
      ylabel('temperature')
      ylim([-0.2 1.2])
      grid on
      M(counter) = getframe(gcf);
      counter = counter + 1;

      % calculate T at subsequent times
      for j=1:tmax
          for i=2:length(T)-1
              T(new,i) = (dt/(dx*dx*p(i)*cp(i)))*(k(i)*(T(old,i+1)-2*T(old,i)+T(old,i-1))...
                                                  + 0.25*(k(i+1)-k(i-1))*(T(old,i+1)-T(old,i-1))) + T(old,i);
          end
          T(old,:) = T(new,:);
          if mod(j,pint-1) == 0
              plot(x,T(old,:),'k','LineWidth',2)
              xlabel('rod')
              ylabel('temperature')
              ylim([-0.2 1.2])
              grid on
              M(counter) = getframe(gcf);
              counter = counter + 1;
          end
      end

      % play and save movie
      f.Visible = 'on';
      v = VideoWriter(sprintf('1-dt=%g',coeff),'MPEG-4');
      v.FrameRate = frate;
      open(v)
      movie(M,plays,frate);
      writeVideo(v,M)
      close(v)
    case 2 % problem 2 (homogeneous only)
      fprintf('Working on Problem 2\n')

      % define temperature array
      Told = zeros(length(x),1);
      Tcur = zeros(length(x),1);

      % ICs
      Told(L/2+1,1) = 1;

      % define thermal conductivity array
      k = zeros(length(x),1);
      k(:) = 1;
      % define density array
      p = zeros(length(x),1);
      p(:) = 1;
      % define specific heat array
      cp = zeros(length(x),1);
      cp(:) = 1;

      % define D s.t. it is as large as possible
      % so that dt is as small as possible
      D = max(k)/(min(p)*min(cp));

      % timestep and tmax
      % coeff > 0.5 does not crash
      coeff = 0.6;
      dt = coeff*dx*dx/D;
      tmax = 150;

      % define s
      s = (D*dt)/(dx*dx);

      % make A matrix
      e = repelem(s,length(x));
      A = spdiags([-e' 1+2*e' -e'],-1:1,length(x),length(x));
      A = full(A);
      
      % initial plotting
      f=figure;
      f.Visible = 'off';
      counter = 1;
      pint = 5;
      plays = 1;
      frate = 4;
      plot(x,Told(:,1)','k','LineWidth',2)
      xlabel('rod')
      ylabel('temperature')
      ylim([-0.2 1.2])
      grid on
      M(counter) = getframe(gcf);
      counter = counter + 1;

      % run simulation tmax times
      for j=1:tmax
          % BCs
          Told(1,1) = 0; Told(length(x),1) = 0;
          Tcur(1,1) = 0; Tcur(length(x),1) = 0;
          % solve linear system for each timestep
          Tcur = A\Told;          
          Told = Tcur;
          % plotting
          if mod(j,pint-1) == 0
              plot(x,Told(:,1),'k','LineWidth',2)
              xlabel('rod')
              ylabel('temperature')
              ylim([-0.2 1.2])
              grid on
              M(counter) = getframe(gcf);
              counter = counter + 1;
          end
      end

      % play and save movie
      f.Visible = 'on';
      v = VideoWriter(sprintf('2-dt=%g',coeff),'MPEG-4');
      v.FrameRate = frate;
      open(v)
      movie(M,plays,frate);
      writeVideo(v,M)
      close(v)
    case 3 % problem 3 (homogeneous only)
      fprintf('Working on Problem 3\n')

      % define temperature array
      Told = zeros(length(x),1);
      Tcur = zeros(length(x),1);

      % ICs
      Told(L/2+1,1) = 1;

      % define thermal conductivity array
      k = zeros(length(x),1);
      k(:) = 1;
      % define density array
      p = zeros(length(x),1);
      p(:) = 1;
      % define specific heat array
      cp = zeros(length(x),1);
      cp(:) = 1;

      % define D s.t. it is as large as possible
      % so that dt is as small as possible
      D = max(k)/(min(p)*min(cp));

      % timestep and tmax
      % coeff > 0.5 does not crash
      coeff = 0.45;
      dt = coeff*dx*dx/D;
      tmax = 150;

      % define s
      s = (D*dt)/(2*dx*dx);

      % make A matrix
      e = repelem(s,length(x));
      A = spdiags([-e' 1+2*e' -e'],-1:1,length(x),length(x));
      A = full(A);

      % make B matrix
      B = spdiags([e' 1-2*e' e'],-1:1,length(x),length(x));
      B = full(B);
      
      % initial plotting
      f=figure;
      f.Visible = 'off';
      counter = 1;
      pint = 5;
      plays = 1;
      frate = 4;
      plot(x,Told(:,1)','k','LineWidth',2)
      xlabel('rod')
      ylabel('temperature')
      ylim([-0.2 1.2])
      grid on
      M(counter) = getframe(gcf);
      counter = counter + 1;

      % run simulation tmax times
      for j=1:tmax
          % BCs
          Told(1,1) = 0; Told(length(x),1) = 0;
          Tcur(1,1) = 0; Tcur(length(x),1) = 0;
          % solve linear system for each timestep
          Tcur = (inv(B)*A)\Told;          
          Told = Tcur;
          % plotting
          if mod(j,pint-1) == 0
              plot(x,Told(:,1),'k','LineWidth',2)
              xlabel('rod')
              ylabel('temperature')
              ylim([-0.2 1.2])
              grid on
              M(counter) = getframe(gcf);
              counter = counter + 1;
          end
      end

      % play and save movie
      f.Visible = 'on';
      v = VideoWriter(sprintf('3-dt=%g',coeff),'MPEG-4');
      v.FrameRate = frate;
      open(v)
      movie(M,plays,frate);
      writeVideo(v,M)
      close(v)
end