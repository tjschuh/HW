function geo441hw05(n,f)
% GEO441HW05(n,f)
%
% simple finite element method (FEM) code for GEO 441 hw 5
%
% INPUT:
%
% n      number of grid points
% f      weight in diffusion eqaution, can be any value
%
% Originally written by tschuh-at-princeton.edu, 03/15/2022
% Last modified by tschuh-at-princeton.edu, 03/17/2022

% number of elements/cells is one less than grid points
ne = n - 1;

% string length
L = 1;

% create grid
x = linspace(0,L,ne+1);

% spacing (doesnt need to be constant)
dx = zeros(1,ne);
for i=1:ne
    dx(i) = x(i+1) - x(i);
end

% BCs
q0 = 1; T1 = 1;

% calculate shape function matrix N
% each row of N corresponds to a different shape function
% each column of N corresponds to position along x array
% N ends up being the identity matrix
N = zeros(ne+1,ne+1);
for i=1:ne+1
    % if on shape function 1
    if i == 1
        N(i,i) = (x(i+1)-x(i))/dx(i);
    % if on shpae function ne+1
    elseif i == ne+1
        N(i,i) = (x(i)-x(i-1))/dx(i-1);
    % if on shape function somewhere in between    
    else
        % there is probably a more elegant way to
        % do this, but analytically we know that
        % it ends up being 1 at one index and 0 elesewhere
        for j=1:ne+1
            if  i == j
                N(i,j) = 1;
            else
                N(i,j) = 0;
            end
        end
    end
end    

% calculate stiffness matrix K
% K has ne x ne dimensions
for i=1:ne
    for j=1:ne
        % if along the diagonal
        if i == j
            % element 1,1 is a special case
            if i ~= 1
                K(i,j) = (1/dx(i-1)) + (1/dx(i));
                % other diagonal elements follow this formula    
            else
                K(i,j) = 1/dx(i);
            end
        % if 1 element away from the diagonal
        elseif abs(i-j) == 1
            K(i,j) = -1/dx(i);
        % all other elements are 0
        else
            K(i,j) = 0;
        end
    end
end

% calculate F
% F has dimensions ne x 1
for i=1:ne
    % 1st element is special because at boundary
    if i == 1
        F(i,1) = q0*N(i,1);
    % middle elements are 0
    elseif i > 1 && i < ne
        F(i,1) = 0;
    % last element is special because at boundary
    else
        F(i,1) = q0*N(i,1) - ((N(i,i+1)-N(i,i))/dx(i))*((N(end,end)-N(end,end-1))/dx(end))*dx(i)*T1;
    end

    % incorporating f into the calculation of F
    % end elements only have 1 shape function so
    % integral is f*area = f*(1/2*dx*h)
    if i == 1 || i == ne
        F(i,1) = F(i,1) + (1/2)*dx(i)*sum(N(i,:))*f;
    % non-end elements have 2 shape functions so
    % integral is f*(1/2*dx*h + 1/2*dx*h)    
    else
        F(i,1) = F(i,1) + (1/2)*dx(i)*sum(N(i,:))*f + (1/2)*dx(i+1)*sum(N(i,:))*f;
    end
end

% calculate d by inverting
d = K\F;

% calculate T
% follow formular given in class
for i=1:ne+1
    T(1,i) = sum(d.*N(1:ne,i)) + T1*N(ne+1,i);
end

% compare to exact solution from class
Tex = T1 + (1 - x).*q0 + (1/2).*(1 - x.^2).*f;

% plotting
p1=plot(x,T,'-b','LineWidth',2);
hold on
scatter(x,T,50,'b','filled')
p2=plot(x,Tex,'--r','LineWidth',2);
scatter(x,Tex,'r','filled')
hold off
longticks([],2)
grid on
xlabel('x')
ylabel('temperature')
title(sprintf('1D steady-state diffusion equation\nnumber of elements = %d, f = %d',n,f))
legend([p1 p2],{'FEM','EXACT'})