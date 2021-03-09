%This is a script to produce a graph for problem 1b

% elastic moduli [GPa]
M1 = 10;
M2 = 100;

x1 = linspace(0,1);
x2 = 1 - x1;

% Voigt bound
MV = x1*M1 + x2*M2;

% Reuss bound
MR = 1./((x1/M1) + (x2/M2));

plot(x1,MV,'b')
xlabel('x_1 volume fraction')
ylabel('M_V [GPa], M_R [GPa]')
ax1=gca;
hold on
plot(x1,MR,'r')
legend('M_V','M_R')
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% set the same Limits and Ticks on ax2 as on ax1;
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
OppTickLabelsX = {'1' '0.9' '0.8' '0.7' '0.6' '0.5' '0.4' '0.3' '0.2' '0.1' '0'};
OppTickLabelsY = {'10' '20' '30' '40' '50' '60' '70' '80' '90' '100'};
% Set the x-tick and y-tick  labels for the second axes
set(ax2, 'XTickLabel', OppTickLabelsX,'YTickLabel',OppTickLabelsY);
xlabel('x_2 volume fraction')
