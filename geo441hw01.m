function geo441hw01()
% GEO441()
%
% code for HW 1
%
% Originally written by tschuh-at-princeton.edu, 02/02/2022
% Last modified by tschuh-at-princeton.edu, 02/06/2022

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

% allocate displacment grid
u = zeros(tmax/dt,xmax/dx);

% Initial Condition
% compute displacement values for first row aka t = 0
for i=1:length(u)
    u(1,i) = exp(-dx*((i/10) - 50)^2);
end

% Dirichlet BCs (fixed ends)
u(:,1) = 0;
u(:,end) = 0;

% compute future times aka subsequent rows by using equation from class
for j=1:tmax/dt-1
    for k=2:xmax/dx-1
        % for row 2 (t = 1), we dont know j - 1 term so skipping that
        if j == 1
            u(j+1,k) = ((c*dt/dx)^2)*(u(j,k+1) - 2*u(j,k) + u(j,k-1)) ...
                       + 2*u(j,k);
        % for every following timestep aka row, we do know j - 1 term    
        else
            u(j+1,k) = ((c*dt/dx)^2)*(u(j,k+1) - 2*u(j,k) + u(j,k-1)) ...
                       + 2*u(j,k) - u(j-1,k);
        end
    end
end

% everything computationally works, but getting wrong values

keyboard