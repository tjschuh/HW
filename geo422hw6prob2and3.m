%saved Frederik's 2 coulumn YEARLY.plt file as a .txt file
%1st column is calendar year, 2nd column is measurement
yearly = importdata('YEARLY.txt');

figure
plot(yearly(:,1),yearly(:,2))
title('YEARLY data')
xlabel('year')
ylabel('mysterious measurement')
xlim([1690 2017])

figure
periodogram(yearly(:,2))
%% 
%saved Frederik's 3 coulumn DECADAL.plt file as a .txt file
%1st column is # of years before calendar year 1950,
%2nd column is measurement, 3rd column is measurement uncertainty
decadal = importdata('DECADAL.txt');

figure
plot(decadal(:,1),decadal(:,2))
title('DECADAL data')
xlabel('# of years before calendar year 1950')
ylabel('mysterious measurement')
xlim([-300 11700])

figure
periodogram(decadal(:,2))