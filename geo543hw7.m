function geo543hw7(N)
% GEO543HW7(N)
%
% Boundary-element code to obtain a numerical solution for the distribution
% of displacements along a strike-slip fault subjected to a uniform stress
% change rather than a uniform displacement discontinuity
%
% INPUT:
%
% N       number of segments to divide fault segment into
%
% OUTPUT:
%
% plot of displacement (s/l) vs position on fault trace (x2/l)
%
% Originally written by tschuh-at-princeton.edu, 10/30/2021

l = 1;
d = l/N;

% set up the correct positions on the fault trace
x = [-l:2*d:l];
x1 = x + d;
x2 = x1(1:N);

% compute the NxN c matrix
for i=1:N                                                             
    for j=1:N                                                             
        c(i,j) = 1/(N*(((2*d*(j-i))/l) - (1/N))*(((2*d*(j-i))/l) + (1/N)));
    end
end

% set the stress elements equal to all ones and we can scale later if needed
T = ones(N,1);

% solve Ax = b
% solve c*s = T
s = linsolve(c,T);

% plotting
scatter(x2,s,'filled','k')
hold on
plot(x2,s,'LineWidth',1.25)
grid on
longticks
xlim([-l l])
xlabel('Position on Fault Trace (x_2/l)')
ylabel('Displacement (s/l)')
