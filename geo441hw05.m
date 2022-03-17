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

% TO-DO:
% fix calculation of N

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

% calculate shape functions
% need to do this properly depsite
% it ending up being identity matrix
N = zeros(ne+1,ne+1);
for i=1:ne+1
    if x(i) < x(2)
        N(1,i) = (x(2)-x(i))/dx(i);
    end
end    

N = eye(ne+1);

% calculate stiffness matrix K
for i=1:ne
    for j=1:ne
        if i == j
            if i ~= 1
                K(i,j) = (1/dx(i-1)) + (1/dx(i));
            else
                K(i,j) = 1/dx(i);
            end
        elseif abs(i-j) == 1
            K(i,j) = -1/dx(i);
        else
            K(i,j) = 0;
        end
    end
end

% calculate F
for i=1:ne
    if i == 1
        F(i,1) = q0*N(i,1);
    elseif i > 1 && i < ne
        F(i,1) = 0;
    else
        F(i,1) = q0*N(i,1) - ((N(i,i+1)-N(i,i))/dx(i))*((N(end,end)-N(end,end-1))/dx(end))*dx(i)*T1;
    end

    if i == 1 || i == ne
        F(i,1) = F(i,1) + (1/2)*dx(i)*sum(N(i,:))*f;
    else
        F(i,1) = F(i,1) + (1/2)*dx(i)*sum(N(i,:))*f + (1/2)*dx(i+1)*sum(N(i,:))*f;
    end
end

% calculate d
d = K\F;

% calculate T
for i=1:ne+1
    T(1,i) = sum(d.*N(1:ne,i)) + T1*N(ne+1,i);
end

% compare to exact solution
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