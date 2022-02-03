function geo441hw01()
% GEO441()
%
% code for HW 1
%
% Originally written by tschuh-at-princeton.edu, 02/02/2022

% wave speed material property for homogenous medium
c = 1;

% grid size
dx = 0.1;

% timestep
dt = dx/c;

% string length
xmax = 100;

% max time
tmax = 10;

%L = [0:dx:xmax];

% allocate displacment grid
u = zeros(tmax/dt,xmax/dx+1);

% initial conditions
for i=1:length(xmax)
    u(1,i) = exp(-dx*(i*dx - (xmax/2))^2);
end
keyboard
%for i=1:xmax/dx

% Dirichlet BCs (fixed ends)
u(1,:) = 0;
u(L,:) = 0;