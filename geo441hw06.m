function geo441hw06()
%
%
% Originally written by tschuh-at-princeton.edu, 03/18/2022
% Last modified by tschuh-at-princeton.edu, 03/21/2022

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
M = zeros(ne,ne);
k = zeros(2,2);
K = zeros(ne,ne);
fe = zeros(2,1);
F = zeros(ne,1);
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
                % this might need to be a + sign
                fe(j,1) = (dx(i)*f/2) - T1*k(j,2);
            else
                fe = (dx(i)*f/2).*[1; 1];
            end
        end
    end
    % now assemble!
    if i ~= ne
        M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + m;
        K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + k;
        F(i:i+1,1) = F(i:i+1,1) + fe;
    else
        M(i,i) = M(i,i) + m(2,2);
        K(i,i) = K(i,i) + k(2,2);
        F(i,1) = F(i,1) + fe(2,1);
    end
end

% now actually perform finite element algorithm
% choose dt and tmax (move to top)
dt = 0.1;
tmax = 10;

% use given initial condition to get d0
d0 = 1 + cos(x);
keyboard
% d0 and K are not same length!
% solve for v0 using d0, M, K, and F
v0 = M\(F-K*d0');

d = zeros(ne,tmax);
v = zeros(ne,tmax);

d(:,1) = d0';
v(:,1) = v0;

dtil = zeros(ne,tmax);
% change i <--> j to be clean
% change last element bc boundary is weird
for j = 1:tmax
    for i = 1:ne
        % prediction of displacement
        dtil(i,j+1) = d(i,j) + (1-alpha)*dt*v(i,j);
    end
    % computation of velocity
    v(:,j+1) = (M + alpha*dt.*K)\(F - K*dtil(:,j+1));
    
    % correction of displacement
    d(:,j+1) = dtil(:,j+1) + alpha*dt*v(:,j+1);
end

% compare against exact solution
for i=1:tmax
    Tex(i,:) = 1 + exp(-i)*cos(x);
end

keyboard