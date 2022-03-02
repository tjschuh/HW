function geo441hw04(n)
% GEO441HW04(n)
%
% Creates and saves a movie of an oscillating 1D string with varying
% boundary conditions and material properties using pseudo-spectral method
%
% INPUT:
%
% n      problem number
%        1 --> problem 1
%        2 --> problem 2
%    
% OUTPUT:
%
% .mp4 file of the movie
%    
% EXAMPLES:
%
% geo441hw04(1)
% geo441hw04(2)
%
% Originally written by tschuh-at-princeton.edu, 02/23/2022
% Last modified by tschuh-at-princeton.edu, 03/02/2022

% grid size and timestep
dx = 0.1;

% string length and max time
L = 100; tmax = 10000;

% create actual string
x = [0:dx:L];

% only ever need 3 rows
% previous timestep --> current timestep --> future timestep
u = zeros(3,L/dx+1);

% allocate velocity grid
v = zeros(3,L/dx+1);

% allocate stress array
T = zeros(3,L/dx+1);

% material properties
switch n
    case 1 % homogeneous
      k = ones(1,L/dx+1); p = ones(1,L/dx+1);
    case 2 % heterogeneous
      split = 60;
      k = ones(1,L/dx+1);
      k(1,(split/dx+1):end) = 4;
      p = ones(1,L/dx+1);
end    
c = sqrt(k./p);

% timestep
%dt = dx/max(c);
dt = 0.01;

% define each row of displacement grid
old = 1; cur = 2; new = 3;

% IC: compute displacement values for t = 0
for j=1:size(u,2)
    u(cur,j) = exp(-0.1*(((j-1)/10) - 50)^2);
end

% displacements for t = -1 are equal to t = 0 values
% need this for 1st timestep
u(old,:) = u(cur,:);

% now compute stress values for t=0 using displacements at t=0
for j=2:size(T,2)-1
    T(cur,j) = (k(j)/(2*dx))*(u(cur,j+1) - u(cur,j-1));
end
T(old,:) = T(cur,:);

% v is always zero at t=0 so we dont need to do any ICs there

% plot t = 0 displacements
f=figure;
f.Visible = 'off';
counter = 1;
pint = 20;
plays = 1;
frate = 12;
plot(x,v(cur,:),'b','LineWidth',2)
hold on
plot(x,T(cur,:),'r','LineWidth',2)
ylim([-1 1])
legend('Velocity','Stress')
grid on
M(counter) = getframe(gcf);
hold off
counter = counter + 1;

grad=gradient(x,dx);
% no propagation, but it works may need to change : --> n element, but
% dont know if I can do fft 1 element at a time
for j=1:tmax
    % compute fft of T and v, multiply by ik
    % shift all elements since fft messes up element order
    % apply inverse fft to get derivative of T and v
    iTf = 2*pi*ifft(grad.*fft(T(cur,:)));
    ivf = 2*pi*ifft(grad.*fft(v(cur,:)));
    
    % use that to find T and v at next timestep
    for m=1:length(T)
        v(new,m) = v(old,m) + (2*dt*iTf(1,m)/p(1,m));
        T(new,m) = T(old,m) + (2*dt*ivf(1,m)*k(1,m));
    end

    % set old = cur, and cur = new
    T(old,:) = T(cur,:); T(cur,:) = T(new,:);
    v(old,:) = v(cur,:); v(cur,:) = v(new,:);

    % make movie and only plot every pint frame
    if mod(j,pint-1) == 0        
        plot(x,real(v(cur,:)),'b','LineWidth',2)
        hold on
        plot(x,real(T(cur,:)),'r','LineWidth',2)
        ylim([-1 1])
        if n == 1
            title('Homogeneous')
        else
            title('Heterogeneous')
        end
        legend('Velocity','Stress')
        grid on
        M(counter) = getframe(gcf);
        hold off
        counter = counter + 1;
    end
end

% play and save movie
f.Visible = 'on';
v = VideoWriter(sprintf('%i',n),'MPEG-4');
v.FrameRate = frate;
open(v)
movie(M,plays,frate);
writeVideo(v,M)
close(v)

function grad=gradient(x,dx)
n=length(x);
val=1/(n*dx);
N=floor((n-1)/2)+1;
p1=[0:N-1];
p2=[-(floor(n/2)):-1];
results=[p2 p1];
grad=i*results*val;
