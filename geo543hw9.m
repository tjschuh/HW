function geo543hw9(N)

% Base Parameters
defval('N',100)
l = 1;
    
% Elliptical crack
d = l/N;
A = zeros(N,N);
B = ones(N,1);
for i=1:N
    for j=1:N
        x2 = 2*d*(j-i);
        c = 1/(((x2/l) + (1/N))*((x2/l) - (1/N)));
        A(i,j) = (1/(N*pi))*c;
    end
end
x = A\B;

x1 = 0:0.1:10;
us_ell = zeros(length(x1),1);

for i=1:length(x1)
    for j=1:N
        x2 = l - ((l/N)*(j - 0.5));
        xplus = (x2./l) + (d/l);
        xminus = (x2./l) - (d/l);
        u = (1/(2*pi))*(atan2(x1(i)./l,xplus) - atan2(x1(i)./l,xminus));
        s_ellipse = -1 + (sum(x(j:N))/N);
        usum(j) = u/s_ellipse;
    end
    us_ell(i,1) = sum(usum,2);
end

% Rectangular crack
s = -1;
d = 1;
x1 = 0:0.1:10;
us_rec = zeros(length(x1),1);

for i=1:length(x1)
    xplus = (x2./l) + (d/l);
    xminus = (x2./l) - (d/l);
    u = (1/(2*pi))*(atan2(x1(i)./l,xplus) - atan2(x1(i)./l,xminus));
    us_rec(i) = u/s;
end

% Scaled Elliptical crack
lprime = l*(4/pi);
d = lprime/N;
A = zeros(N,N);
B = ones(N,1);
for i=1:N
    for j=1:N
        x2 = 2*d*(j-i);
        c = 1/(((x2/lprime) + (1/N))*((x2/lprime) - (1/N)));
        A(i,j) = (1/(N*pi))*c;
    end
end
x = A\B;

x1 = 0:0.1:10;
us_sell = zeros(length(x1),1);

for i=1:length(x1)
    for j = 1:N
        x2 = lprime - ((lprime/N)*(j - 0.5));
        xplus = (x2./lprime) + (d/lprime);
        xminus = (x2./lprime) - (d/lprime);
        u = (1/(2*pi))*(atan2(x1(i)./lprime,xplus) - atan2(x1(i)./lprime,xminus));
        s_ellipse = -1 + (sum(x(j:N))/N);
        usum(j) = u/s_ellipse;
    end
    us_sell(i,1) = sum(usum,2);
end

% Plotting
plot(x1./l,us_ell,'r','LineWidth',1.25)
hold on
plot(x1./l,us_rec,'g','LineWidth',1.25)
plot(x1./l,us_sell,'b','LineWidth',1.25)
legend('Ellispe','Rectangle','Scaled Ellispe');
longticks
grid on
xlabel('Distance from Fault Trace (x_1/l)')
ylabel('Surface Displacement (\mu/s)')