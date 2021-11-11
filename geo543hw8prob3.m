function geo543hw8prob3()

% calculate slip for the ellipse 
N = 30;
Nmid = round(N/2);
slip = zeros(N,1);
% Tau
b = ones(N,1);
A = zeros(N,N);
% hold the x2 values
x2 = zeros(N,N);
% defined from -1 to +1 
l2 = 2;
l = 1; 
d = l/N;
% midpoint of each element
d2 = l2/N; 
for j = 1:N
    for i = 1:N
        x2(j,i) = d2*(j-i);
        A(j,i) = (1/(pi*N))/((x2(j,i)/l-1/N)*(x2(j,i)/l+1/N));
    end
end
slip = linsolve(A,b);
dist = x2(Nmid,:);

% vary x2 for different values of d 
U3acc_e = 0;
U3acc_r = 0;
for j = 1:N*100
    x1 = (j-1)/(100*l);
    for i = 1:Nmid
        % only for odd numbers
        x2 = d*i;
        U3temp_e = (slip(Nmid+i)/(2*pi))*(atan2(x1,dist(Nmid)/d)-atan2(x1,dist(Nmid)-d))/(-sum(slip));
        U3acc_e = U3acc_e + U3temp_e;
        U3temp_r = (mean(slip(Nmid+i,:))/(2*pi))*(atan2(x1,d)-atan2(x1,-d))/(-sum(slip));
        U3acc_r = U3acc_r + U3temp_r;
    end
    %only using half the space 
    U3_ellispe(j) = U3acc_e*2;
    U3_regular(j) = U3acc_r*2;
    U3acc_e = 0;
    U3acc_r = 0;
    U3xe(j) = x1;
    U3xes(j) = x1*(4/pi);
    U3xr(j) = x1;
end

% plotting
plot(U3xe,U3_ellispe,'LineWidth',1.25);
hold on
plot(U3xes,U3_ellispe,'LineWidth',1.25);
plot(U3xr,U3_regular,'LineWidth',1.25);
legend('Scaled Ellispe','Regular','Ellispe');
xlabel('Distance from Fault')
ylabel('Surface Displacement')
xlim([0 0.6])
ylim([0 0.6])
longticks