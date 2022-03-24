function geo441hw06(prob,n)
% GEO441HW06(prob,n)
%
% Creates and saves a movie file of Finite Element Method
% (FEM) solution of 1d unsteady-state diffusion equation
%
% INPUT:
%
% prob     which problem you want to solve [a or b]
% n        number of grid points in your FEM scheme
%
% OUTPUT:
%
% .mp4 movie file
%
% EXAMPLE:
%
% geo441hw06('a',10)
%
% Originally written by tschuh-at-princeton.edu, 03/18/2022
% Last modified by tschuh-at-princeton.edu, 03/23/2022

% number of elements/cells
ne = n - 1;

switch prob
    case 'a'
      fprintf('Solving Problem a...\n')
      % spatial length
      L = pi/2;
      % material properties
      p = 1; cp = 1; kap = 1;
      % BCs
      q0 = 0; TL = 1;
      % force term
      f = 0;
      % timestep and simulation length
      dt = 0.01; tmax = 200;
      % Crank-Nicolson scheme
      alpha = 0.5;
    case 'b'
      fprintf('Solving Problem b...\n')
      L = 20;
      p = 1; cp = 1; kap = 1;
      q0 = 0; TL = 0; Tm = 1;
      f = 0;
      dt = 0.01; tmax = 200;
      alpha = 0.5;
    otherwise
      error('Please enter a valid problem (a or b)')
end

% create grid
x = linspace(0,L,ne+1);

% spacing (doesnt need to be constant)
dx = zeros(1,ne);
for i=1:ne
    dx(i) = x(i+1) - x(i);
end

% local to global scheme
% use local 2x2 N, m, and k matrices to build
% global (ne+1)x(ne+1) M and K matrices where
% N=shape functions, M=capacity matrix, and K=stiffness matrix
% use local 2x1 fe RHS vector to build global (ne+1)x1 F RHS vector
N = zeros(2,2);
m = zeros(2,2);
M = zeros(ne+1,ne+1);
k = zeros(2,2);
K = zeros(ne+1,ne+1);
fe = zeros(2,1);
F = zeros(ne+1,1);
% i=globael domain
% j=local domain
for i=1:ne
    % define xA and xA+1 from global coordinates
    xA = x(i); xA1 = x(i+1);
    for j=1:2
        % b=column number of local matrices
        for b=1:2
            if j == 1
                xp = xA;
            else
                xp = xA1;
            end
            % choose linear mapping between global and local
            e(j) = (2*xp-xA-xA1)/dx(i);
            % calculate N
            % we want 2x2 identity matrix essentially
            if b == 1
                N(j,b) = 0.5*(1+e(j)*-1);
            else
                N(j,b) = 0.5*(1+e(j)*1);
            end
            % calculate derivatives of N
            dN(j) = e(j)/2;
            % calculate m and k
            % m = (p*cp*dx/6).*[2 1; 1 2]
            % k = (kap/dx).*[1 -1; -1 1]
            if j == b
                m(j,b) = (dx(i)*p*cp/6)*2;
                k(j,b) = (kap/dx(i))*(2*1/2);
            else
                m(j,b) = (dx(i)*p*cp/6)*1;
                k(j,b) = (kap/dx(i))*(2*-1/2);
            end
            % calculate f using BCs, force term, N, and k
            if i == 1
                fe(j,1) = (dx(i)*f/2) + q0*N(j,1);
            elseif i == ne
                fe(j,1) = (dx(i)*f/2) - TL*k(j,2);
            else
                fe = (dx(i)*f/2).*[1; 1];
            end
        end
    end
    % now assemble by placing 2x2 matrices
    % into global (ne+1)x(ne+1) matrices
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + m;
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + k;
    % last element of F will be nonzero, but we only care
    % about elements 1:end-1 since F has elements ne+1
    F(i:i+1,1) = F(i:i+1,1) + fe;
end

% now actually perform finite element algorithm (FEM)
% use given initial condition to get d0
% d = T in this code
switch prob
    case 'a'
      % start at 2, end at 1
      d0 = 1 + cos(x);
    case 'b'
      % [1 1 1 ... 1 1 0]
      d0 = [ones(1,length(x)-1) TL];
end

% estimate v0 using d0, M, K, and F and inverting
% but only use elements 1:end-1 b/c we want to keep T=1 at boundary
v0 = (M(1:end-1,1:end-1)+alpha*dt.*K(1:end-1,1:end-1))\(F(1:end-1,1)-K(1:end-1,1:end-1)*d0(1,1:end-1)');

% allocate d, dtil, and v arrays
d = zeros(2,ne+1);
dtil = zeros(1,ne+1);
v = zeros(2,ne+1);

% for convenience
old = 1; new = 2;

% first entries of d and v equal d0 and v0
d(old,:) = d0;
v(old,1:end-1) = v0';

switch prob
    case 'a'
      Tex(1,:) = 1 + cos(x);
    case 'b'
      % erfc(t=0) = erfc(inf) = 0
      % Tex(1,:) = [1 1 1 ... 1 1 0]
      Tex(1,:) = Tm*[1-zeros(1,length(x)-1) TL];
end
      
% initial plotting
f=figure;
f.Visible = 'off';
counter = 1;
% plotting interval
pint = 5;
% number of movie plays
plays = 1;
% frame rate
frate = 4;
% scatter plot size
sz = 50;
p1=plot(x,d(old,:),'-b','LineWidth',2);
hold on
scatter(x,d(old,:),sz,'b','filled')
p2=plot(x,Tex(1,:),'--r','LineWidth',2);
scatter(x,Tex(1,:),sz,'r','filled')
title(sprintf('1D unsteady-state diffusion equation\nnumber of grid points = %d, problem %s',n,prob))
xlabel('x')
ylabel('temperature')
switch prob
    case 'a'
      ylim([1 2])
    case 'b'
      ylim([0 1.2])
end
legend([p1 p2],{'FEM','EXACT'})
grid on
% save movie frame
GF(counter) = getframe(gcf);
hold off
counter = counter + 1;

% predictor-corrector integration scheme
for i=1:tmax
    % prediction of displacement
    dtil(1,:) = d(old,:) + (1-alpha)*dt*v(old,:);

    % computation of velocity
    % ignoring last element when inverting b/c T=1 at boundary
    v(new,1:end-1) = (M(1:end-1,1:end-1) + alpha*dt.*K(1:end-1,1:end-1))\(F(1:end-1,1) - K(1:end-1,1:end-1)*dtil(1,1:end-1)');
    
    % correction of displacement
    d(new,:) = dtil(1,:) + alpha*dt.*v(new,:);

    % compare against exact solution
    switch prob
        case 'a'
          Tex(1,:) = 1 + exp(-i*dt)*cos(x);
        case 'b'
          %Tex(1,:) = Tm*(1 - erfc(x./(2*sqrt((kap*i*dt)/(p*cp)))));
          Tex(1,:) = 1 - erfc((L-x)./(2*sqrt((kap*i*dt)/(p*cp))));
    end

    % set row 1 (old) = row 2 (new) and plot row 1 (old)
    d(old,:) = d(new,:);
    v(old,:) = v(new,:);

    % plotting
    if mod(i,pint-1) == 0
        p1=plot(x,d(old,:),'-b','LineWidth',2);
        hold on
        scatter(x,d(old,:),sz,'b','filled')
        p2=plot(x,Tex(1,:),'--r','LineWidth',2);
        scatter(x,Tex(1,:),sz,'r','filled')
        title(sprintf('1D unsteady-state diffusion equation\nnumber of grid points = %d, problem %s',n,prob))
        xlabel('x')
        ylabel('temperature')
        %ylim([1 2])
        switch prob
          case 'a'
            ylim([1 2])
          case 'b'
            ylim([0 1.2])
        end
        legend([p1 p2],{'FEM','EXACT'})
        grid on
        GF(counter) = getframe(gcf);
        hold off
        counter = counter + 1;
    end
end

% play and save movie
f.Visible = 'on';
v = VideoWriter(sprintf('problem-%s',prob),'MPEG-4');
v.FrameRate = frate;
open(v)
movie(GF,plays,frate);
writeVideo(v,GF)
close(v)
fprintf('DONE\n')