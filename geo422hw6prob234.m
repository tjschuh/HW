%Problem 2
%saved Frederik's 2 coulumn YEARLY.plt file as a .txt file
%1st column is calendar year, 2nd column is measurement
yearly = importdata('YEARLY.txt');

subplot(2,3,1)
plot(yearly(:,1),yearly(:,2))
title('YEARLY data')
xlabel('year')
ylabel('mysterious measurement')
xlim([1690 2017])

subplot(2,3,2)
periodogram(yearly(:,2))
%% 
%Problem 3
%saved Frederik's 3 coulumn DECADAL.plt file as a .txt file
%1st column is # of years before calendar year 1950,
%2nd column is measurement, 3rd column is measurement uncertainty
decadal = importdata('DECADAL.txt');

subplot(2,3,4)
plot(decadal(:,1),decadal(:,2))
title('DECADAL data')
xlabel('# of years before calendar year 1950')
ylabel('mysterious measurement')
xlim([-300 11700])

subplot(2,3,5)
periodogram(decadal(:,2))
%% 
%Problem 4
%popular data windows are:
%bartlett, blackman, chebwin, hamming, hann, kaiser

subplot(2,3,3)
periodogram(yearly(:,2),hann(length(yearly(:,2))))

subplot(236)
periodogram(decadal(:,2),kaiser(length(decadal(:,2))))
