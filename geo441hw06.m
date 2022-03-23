function geo441hw06()
%
%
% Originally written by tschuh-at-princeton.edu, 03/18/2022
% Last modified by tschuh-at-princeton.edu, 03/23/2022

% number of grid points
n = 10;

% number of elements/cells
ne = n - 1;

% string length
L = pi/2;

% create grid
x = linspace(0,L,ne+1);

% spacing (doesnt need to be constant)
dx = zeros(1,ne);
for i=1:ne
    dx(i) = x(i+1) - x(i);
end

% timestep and simulation length
dt = 0.01; tmax = 200;

% material properties
p = 1; cp = 1; kap = 1;

% BCs
q0 = 0; T1 = 1;

% force term
f = 0;

% Crank-Nicolson scheme
alpha = 0.5;

N = zeros(2,2);
m = zeros(2,2);
M = zeros(ne+1,ne+1);
k = zeros(2,2);
K = zeros(ne+1,ne+1);
fe = zeros(2,1);
F = zeros(ne+1,1);
for i=1:ne
    xA = x(i); xA1 = x(i+1);
    for j=1:2
        for b=1:2
            if j == 1
                xp = xA;
            else
                xp = xA1;
            end
            e(j) = (2*xp-xA-xA1)/dx(i);
            % calculate N
            if b == 1
                N(j,b) = 0.5*(1+e(j)*-1);
            else
                N(j,b) = 0.5*(1+e(j)*1);
            end
            % calculate derivatives
            dN(j) = e(j)/2;
            % calculate m and k
            if j == b
                m(j,b) = (dx(i)*p*cp/6)*2;
                k(j,b) = (kap/dx(i))*(2*1/2);
            else
                m(j,b) = (dx(i)*p*cp/6)*1;
                k(j,b) = (kap/dx(i))*(2*-1/2);
            end
            % calculate f now
            if i == 1
                fe(j,1) = (dx(i)*f/2) + q0*N(j,1);
            elseif i == ne
                fe(j,1) = (dx(i)*f/2) - T1*k(j,2);
            else
                fe = (dx(i)*f/2).*[1; 1];
            end
        end
    end
    % now assemble!
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + m;
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + k;
    % last element of F will be nonzero, but we only care
    % about elements 1:end-1 since F has elements ne+1
    F(i:i+1,1) = F(i:i+1,1) + fe;
end

% now actually perform finite element algorithm
% use given initial condition to get d0
% d = T
d0 = 1 + cos(x);

% solve for v0 using d0, M, K, and F
% but only use elements up to last bc we want
% to keep T=1 at boundary
v0 = (M(1:end-1,1:end-1)+alpha*dt.*K(1:end-1,1:end-1))\(F(1:end-1,1)-K(1:end-1,1:end-1)*d0(1,1:end-1)');

d = zeros(tmax,ne+1);
v = zeros(tmax,ne+1);

d = zeros(2,ne+1);
dtil = zeros(1,ne+1);
v = zeros(2,ne+1);

old = 1; new = 2;

d(old,:) = d0;
v(old,1:end-1) = v0';

Tex(1,:) = 1 + cos(x);

% initial plotting
% add title with n and scatter plot
f=figure;
f.Visible = 'off';
counter = 1;
pint = 5;
plays = 1;
frate = 4;
sz = 50;
p1=plot(x,d(old,:),'-b','LineWidth',2);
hold on
scatter(x,d(old,:),sz,'b','filled')
p2=plot(x,Tex(1,:),'--r','LineWidth',2);
scatter(x,Tex(1,:),sz,'r','filled')
title(sprintf('1D unsteady-state diffusion equation\nnumber of grid points = %d',n))
xlabel('x')
ylabel('temperature')
ylim([1 2])
legend([p1 p2],{'FEM','EXACT'})
grid on
GF(counter) = getframe(gcf);
hold off
counter = counter + 1;

% ignored last element bc T=1 at boundary
for i=1:tmax
    % prediction of displacement
    dtil(1,:) = d(old,:) + (1-alpha)*dt*v(old,:);

    % computation of velocity
    v(new,1:end-1) = (M(1:end-1,1:end-1) + alpha*dt.*K(1:end-1,1:end-1))\(F(1:end-1,1) - K(1:end-1,1:end-1)*dtil(1,1:end-1)');
    
    % correction of displacement
    d(new,:) = dtil(1,:) + alpha*dt.*v(new,:);

    % compare against exact solution
    Tex(1,:) = 1 + exp(-i*dt)*cos(x);

    d(old,:) = d(new,:);
    v(old,:) = v(new,:);

    % plotting
    if mod(i,pint-1) == 0
        p1=plot(x,d(old,:),'-b','LineWidth',2);
        hold on
        scatter(x,d(old,:),sz,'b','filled')
        p2=plot(x,Tex(1,:),'--r','LineWidth',2);
        scatter(x,Tex(1,:),sz,'r','filled')
        title(sprintf('1D unsteady-state diffusion equation\nnumber of grid points = %d',n))
        xlabel('x')
        ylabel('temperature')
        ylim([1 2])
        legend([p1 p2],{'FEM','EXACT'})
        grid on
        GF(counter) = getframe(gcf);
        hold off
        counter = counter + 1;
    end
end

% play and save movie
f.Visible = 'on';
v = VideoWriter(sprintf('a'),'MPEG-4');
v.FrameRate = frate;
open(v)
movie(GF,plays,frate);
writeVideo(v,GF)
close(v)