function geo441hw04()
% GEO441HW04()
%
% Creates and saves a movie of an oscillating 1D string with varying
% boundary conditions and material properties using pseudo-spectral method
%
% INPUT:
%    
% OUTPUT:
%
% .mp4 file of the movie
%    
% Originally written by tschuh-at-princeton.edu, 02/23/2022

% material properties
% wave speed c = sqrt(k/p)
k = 1; p = 1;
c = sqrt(k/p);

% grid size and timestep
dx = 0.1; dt = dx/c;

% string length and max time
L = 100; tmax = 2000;

% create actual string
x = [0:dx:L];

% only ever need 3 rows
% previous timestep --> current timestep --> future timestep
u = zeros(3,L/dx+1);

% allocate velocity grid
v = zeros(3,L/dx+1);

% allocate stress array
T = zeros(3,L/dx+1);

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

keyboard

% compute fft of T, multiply by ik, and shift all the elements since fft messes up element order
Tf = i*((2*pi)/(dx*length(x)))*fft(T(cur,:));

% apply inverse fft to get derivative of T
% use that to find v at next timestep
% do same thing with v to find next T