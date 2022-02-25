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

% grid size and timestep
dx = 0.1;

% string length and max time
L = 100; tmax = 300;

% create actual string
x = [0:dx:L];

% only ever need 3 rows
% previous timestep --> current timestep --> future timestep
u = zeros(3,L/dx+1);

% allocate velocity grid
v = zeros(3,L/dx+1);

% allocate stress array
T = zeros(3,L/dx+1);

% material properties (homogeneous for now)
k = ones(1,L/dx+1); p = ones(1,L/dx+1);
c = sqrt(k./p);

% timestep
dt = dx/max(c);
dt = 0.001;

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
pint = 5;
plays = 1;
frate = 15;      
plot(x,v(cur,:),'b','LineWidth',2)
hold on
plot(x,T(cur,:),'r','LineWidth',2)
ylim([-1 1])
legend('Velocity','Stress')
grid on
M(counter) = getframe(gcf);
hold off
counter = counter + 1;

multi=grad(x,dx);
% no propagation, but it works
for j=1:tmax
    % compute fft of T and v, multiply by ik
    % shift all elements since fft messes up element order
    Tf = multi.*fft(T(cur,:));
    vf = multi.*fft(v(cur,:));

    % apply inverse fft to get derivative of T and v
    iTf = ifft(Tf);
    ivf = ifft(vf);

    % use that to find T and v at next timestep
    T(new,:) = T(old,:) + (2*dt.*k(1,:)).*ivf(1,:);
    v(new,:) = v(old,:) + (2*dt./p(1,:)).*iTf(1,:);
    
    % set old = cur, and cur = new
    T(old,:) = T(cur,:); T(cur,:) = T(new,:);
    v(old,:) = v(cur,:); v(cur,:) = v(new,:);

    % make movie and only plot every pint frame
    if mod(j,pint-1) == 0        
        plot(x,real(v(cur,:)),'b','LineWidth',2)
        hold on
        plot(x,real(T(cur,:)),'r','LineWidth',2)
        ylim([-1 1])
        title('Homogeneous')
        legend('Velocity','Stress')
        grid on
        M(counter) = getframe(gcf);
        hold off
        counter = counter + 1;
    end
end

% play and save movie
f.Visible = 'on';
v = VideoWriter('1','MPEG-4');
v.FrameRate = frate;
open(v)
movie(M,plays,frate);
writeVideo(v,M)
close(v)

function multi=grad(x,dx)
n=length(x);
val=(2*pi)/(n*dx);
N=floor((n-1)/2)+1;
p1=[0:N-1];
p2=[-(floor(n/2)):-1];
results=[p2 p1];
blu=i*results*val;
