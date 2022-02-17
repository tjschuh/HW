function geo441hw03()
%
%
% Originally written by tschuh-at-princeton.edu, 02/17/2022

% grid size and rod length
dx = 1; L = 100;

% create actual rod, N = (L/dx) + 1
x = [0:dx:L];

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