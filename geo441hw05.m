function geo441hw05()
%
% REQUIRES:
%
% Symbolic Math Toolbox
%
% Originally written by tschuh-at-princeton.edu, 03/15/2022

% number of elements
n = 3;

% string length
L = 1;

% spatial scale
dx = 1/(n-1);

% create grid
x = 0:dx:L;

% for i=1:length(x)
%     if x(i) < dx
%         N1(i) = 1 - 2*x(i);
%         N2(i) = 2*x(i);
%         N3(i) = 0;
%     else
%         N1(i) = 0;
%         N2(i) = 2*(1 - x(i));
%         N3(i) = 2*x(i) - 1;
%     end
% end

syms x
N1A = 1 - 2*x;
N1B = 0;
N2A = 2*x;
N2B = 2*(1 - x);
N3A = 0;
N3B = 2*x - 1;

N = [N1A N1B; N2A N2B; N3A N3B];

for j=1:n-1
    for k=1:n-1
        K(j,k) = (1/2)*diff(N(j,1))*diff(N(k,1)) + (1/2)*diff(N(j,2))*diff(N(k,2));
    end
end

keyboard