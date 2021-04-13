% This is a script for GEO 424 HW 9

% Originally written by Terance Schuh, 4/12/2021

% load in ISC data that is given
load ISC_data

% distance in degrees up to 24
deltad = ISC_data(:,1);

% convert delta to radians
deltar = (pi/180)*deltad; 

% travel times in seconds
T = ISC_data(:,2);

% computer ray parameter p, p = dT/d(delta)
for i = 2:size(T,1)
    p(i-1) = (T(i)-T(i-1))/(deltar(i)-deltar(i-1));
end
p = p';

figure
plot(deltar,T)
xlabel('distance{\it \Delta} [radians]')
ylabel('travel time{\it T} [s]')
grid on

figure
plot(deltar(2:end),p)
xlabel('distance{\it \Delta} [radians]')
ylabel('ray parameter{\it p}')
grid on

% radius of earth [km]
r0 = 6371;
dx = 2;
eta = p;
counter = 1;

for L = 0:2:96
  if L == 0
    r(counter,1) = r0;
    counter = counter + 1;
  else
    N = (L/dx) + 1;
    for j = 2:N
      f1 = log(p(j-1,1)/eta(N,1) + ((p(j-1,1)/eta(N,1))^2 - 1)^(1/2));
      f2 = log(p(j,1)/eta(N,1) + ((p(j,1)/eta(N,1))^2 - 1)^(1/2));
      trap(j-1) = ((f1+f2)*dx)/2;
    end
    area = sum(trap);
    r(counter,1) = r0*exp((-1/pi)*area);
    counter = counter + 1;
  end
end

% velocity v [km/s] --> r [km], eta [s]
v = r./eta;

figure
plot(r,v)
grid on
xlabel('radius{\it r} [km]')
ylabel('velocity{\it v} [km/s]')
