% define lambda's
modlambda238=1.55125e-10;
modlambda235=9.8485e-10;

% create matrix of data from problem set 
UPbt0238=[18;14;100;7.5];
Pbt0206=9.307*ones(4,1);
Pbt0207=10.294*ones(4,1);
Data=[UPbt0238 Pbt0206 Pbt0207];
Times=[0, 0.1e9, 1e9, 3e9, 4.567e9];


%Calculate 238/235 U at each time point
%First, find the initial 238U using modern ratio
U238o=137.82/(exp(-modlambda238*4.567e9));
Pb204=U238o./Data(:,1);

%then find 238U at each time point
for index=1:1:5
    t=Times(1,index);
    U238(1,index)=U238o*exp(-modlambda238*t);
end
%Find 238U/204Pb ratios at each time point 
U238Pb204Matrix=U238./Pb204; %this is a matrix each row different compenent


%do the same for 235U
U235o=1/(exp(-modlambda235*4.567e9));
for index=1:1:5
    t=Times(1,index);
    U235(1,index)=U235o*exp(-modlambda235*t);
end
U235Pb204Matrix=U235./Pb204;

%Calculate 206/204 values at each time for each material
for index=1:1:5
    t=Times(1,index);
    Pb206(:,index)=Pbt0206+(U238Pb204Matrix(:,index)*(exp(modlambda238*t)-1));
end

%Calculate 207/204 values at each time for each material 
for index=1:1:5
    t=Times(1,index);
    Pb207(:,index)=Pbt0207+(U235Pb204Matrix(:,index)*(exp(modlambda235*t)-1));
end

%Calculate values that will be plotted
Pb204206=Pb206.^-1;
Pb207206=Pb207.*Pb204206;

figure
plot(Pb204206(:,1),Pb207206(:,1),'-xc'); %t=0
hold on 
plot(Pb204206(:,2),Pb207206(:,2),'-xr'); %t=100 Ma
plot(Pb204206(:,3),Pb207206(:,3),'-xb'); %t=1 Ga
plot(Pb204206(:,4),Pb207206(:,4),'-xg'); %t=3 Ga
plot(Pb204206(:,5),Pb207206(:,5),'-xk'); %t=4.567 Ga
legend('t=100 Ma','t=1 Ga', 't=3 Ga', 't=4.567 Ga','Location', 'northwest')
xlabel('204 Pb/206 Pb','FontSize', 14)
ylabel('207 Pb/206 Pb','FontSize', 14)